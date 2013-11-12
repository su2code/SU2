/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.8
 *
 * Stanford University Unstructured (SU2).
 * Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
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


COutput::COutput(void) {
  
  /*--- Initialize point and connectivity counters to zero. ---*/
  nGlobal_Poin      = 0;
  nSurf_Poin        = 0;
  nGlobal_Elem      = 0;
  nSurf_Elem        = 0;
  nGlobal_Tria      = 0;
  nGlobal_Quad      = 0;
  nGlobal_Tetr      = 0;
  nGlobal_Hexa      = 0;
  nGlobal_Wedg      = 0;
  nGlobal_Pyra      = 0;
  nGlobal_Line      = 0;
  nGlobal_BoundTria = 0;
  nGlobal_BoundQuad = 0;
  
  /*--- Initialize CGNS write flag ---*/
  wrote_base_file = false;
  
  /*--- Initialize CGNS write flag ---*/
  wrote_CGNS_base = false;
  
  /*--- Initialize Tecplot write flag ---*/
  wrote_Tecplot_base = false;
  
  /*--- Initialize Paraview write flag ---*/
  wrote_Paraview_base = false;
  
}

COutput::~COutput(void) { }

void COutput::SetSurfaceCSV_Flow(CConfig *config, CGeometry *geometry,
                                 CSolver *FlowSolver, unsigned long iExtIter,
                                 unsigned short val_iZone) {
  
  unsigned short iMarker;
  unsigned long iPoint, iVertex, Global_Index;
  double PressCoeff = 0.0, SkinFrictionCoeff, HeatFlux, xCoord, yCoord, zCoord, Mach, Pressure;
  char cstr[200];
  
  unsigned short solver = config->GetKind_Solver();
  unsigned short nDim = geometry->GetnDim();
  
#ifdef NO_MPI
  
  char buffer [50];
  ofstream SurfFlow_file;
  
  
  /*--- Write file name with extension if unsteady ---*/
  strcpy (cstr, config->GetSurfFlowCoeff_FileName().c_str());
  
  if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.csv", int(val_iZone));
    if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.csv", int(val_iZone));
    if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.csv", int(val_iZone));
    if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.csv", int(val_iZone));
    if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.csv", int(val_iZone));
    
  } else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.csv", int(iExtIter));
    if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.csv",  int(iExtIter));
    if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.csv",   int(iExtIter));
    if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv",    int(iExtIter));
    if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
  }
  else
    sprintf (buffer, ".csv");
  
  strcat (cstr, buffer);
  SurfFlow_file.precision(15);
  SurfFlow_file.open(cstr, ios::out);
  
  SurfFlow_file << "\"Global_Index\", \"x_coord\", \"y_coord\", ";
  if (nDim == 3) SurfFlow_file << "\"z_coord\", ";
  SurfFlow_file << "\"Pressure\", \"Pressure_Coefficient\", ";

  switch (solver) {
    case EULER : SurfFlow_file <<  "\"Mach_Number\"" << endl; break;
    case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"Skin_Friction_Coefficient\"" << endl; break;
    case TNE2_EULER: SurfFlow_file << "\"Mach_Number\"" << endl; break;
    case TNE2_NAVIER_STOKES: SurfFlow_file << "\"Skin_Friction_Coefficient\", \"Heat_Flux\"" << endl; break;
  }
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        xCoord = geometry->node[iPoint]->GetCoord(0);
        yCoord = geometry->node[iPoint]->GetCoord(1);
        if (nDim == 3) zCoord = geometry->node[iPoint]->GetCoord(2);
        Pressure = FlowSolver->node[iPoint]->GetPressure(COMPRESSIBLE);
        PressCoeff = FlowSolver->GetCPressure(iMarker,iVertex);
        SurfFlow_file << scientific << Global_Index << ", " << xCoord << ", " << yCoord << ", ";
        if (nDim == 3) SurfFlow_file << scientific << zCoord << ", ";
        SurfFlow_file << scientific << Pressure << ", " << PressCoeff << ", ";
        switch (solver) {
          case EULER :
            Mach = sqrt(FlowSolver->node[iPoint]->GetVelocity2()) / FlowSolver->node[iPoint]->GetSoundSpeed();
            SurfFlow_file << scientific << Mach << endl;
            break;
          case NAVIER_STOKES: case RANS:
            SkinFrictionCoeff = FlowSolver->GetCSkinFriction(iMarker,iVertex);
            SurfFlow_file << scientific << SkinFrictionCoeff << endl;
            break;
          case TNE2_EULER:
            Mach = sqrt(FlowSolver->node[iPoint]->GetVelocity2()) / FlowSolver->node[iPoint]->GetSoundSpeed();
            SurfFlow_file << scientific << Mach << endl;
            break;
          case TNE2_NAVIER_STOKES:
            SkinFrictionCoeff = FlowSolver->GetCSkinFriction(iMarker,iVertex);
            HeatFlux = FlowSolver->GetHeatTransferCoeff(iMarker,iVertex);
            SurfFlow_file << scientific << SkinFrictionCoeff << ", " << HeatFlux << endl;
        }
      }
    }
  }
  
  SurfFlow_file.close();
  
#else
  
  int rank = MPI::COMM_WORLD.Get_rank(), iProcessor, nProcessor = MPI::COMM_WORLD.Get_size();

  unsigned long Buffer_Send_nVertex[1], nVertex_Surface = 0, nLocalVertex_Surface = 0,
  MaxLocalVertex_Surface = 0, nBuffer_Scalar, *Buffer_Receive_nVertex = NULL, position;
  ofstream SurfFlow_file;
  
  /*--- Write the surface .csv file, the information of all the vertices is
   sent to the MASTER_NODE for writing ---*/
  nLocalVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface++;
      }
  
  if (rank == MASTER_NODE)
    Buffer_Receive_nVertex = new unsigned long [nProcessor];
  
  Buffer_Send_nVertex[0] = nLocalVertex_Surface;
  
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI::UNSIGNED_LONG, MPI::MAX);
  MPI::COMM_WORLD.Gather(&Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
  
  double *Buffer_Send_Coord_x = new double [MaxLocalVertex_Surface];
  double *Buffer_Send_Coord_y = new double [MaxLocalVertex_Surface];
  double *Buffer_Send_Coord_z = new double [MaxLocalVertex_Surface];
  double *Buffer_Send_Press = new double [MaxLocalVertex_Surface];
  double *Buffer_Send_CPress = new double [MaxLocalVertex_Surface];
  double *Buffer_Send_Mach = new double [MaxLocalVertex_Surface];
  double *Buffer_Send_SkinFriction = new double [MaxLocalVertex_Surface];
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long [MaxLocalVertex_Surface];
  
  nVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          Buffer_Send_Press[nVertex_Surface] = FlowSolver->node[iPoint]->GetPressure(COMPRESSIBLE);
          Buffer_Send_CPress[nVertex_Surface] = FlowSolver->GetCPressure(iMarker,iVertex);
          Buffer_Send_Coord_x[nVertex_Surface] = geometry->node[iPoint]->GetCoord(0);
          Buffer_Send_Coord_y[nVertex_Surface] = geometry->node[iPoint]->GetCoord(1);
          Buffer_Send_GlobalIndex[nVertex_Surface] = geometry->node[iPoint]->GetGlobalIndex();
          if (nDim == 3) Buffer_Send_Coord_z[nVertex_Surface] = geometry->node[iPoint]->GetCoord(2);
          if (config->GetKind_Solver() == EULER)
            Buffer_Send_Mach[nVertex_Surface] = sqrt(FlowSolver->node[iPoint]->GetVelocity2()) / FlowSolver->node[iPoint]->GetSoundSpeed();
          if ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver() == RANS))
            Buffer_Send_SkinFriction[nVertex_Surface] = FlowSolver->GetCSkinFriction(iMarker,iVertex);
          nVertex_Surface++;
        }
      }
  
  double *Buffer_Receive_Coord_x = NULL, *Buffer_Receive_Coord_y = NULL, *Buffer_Receive_Coord_z = NULL, *Buffer_Receive_Press = NULL, *Buffer_Receive_CPress = NULL,
  *Buffer_Receive_Mach = NULL, *Buffer_Receive_SkinFriction = NULL;
  unsigned long *Buffer_Receive_GlobalIndex = NULL;
  
  if (rank == MASTER_NODE) {
    Buffer_Receive_Coord_x = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Coord_y = new double [nProcessor*MaxLocalVertex_Surface];
    if (nDim == 3) Buffer_Receive_Coord_z = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Press = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_CPress = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Mach = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_SkinFriction = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_GlobalIndex = new  unsigned long [nProcessor*MaxLocalVertex_Surface];
  }
  
  nBuffer_Scalar = MaxLocalVertex_Surface;
  
  /*--- Send the information to the Master node ---*/
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Gather(Buffer_Send_Coord_x, nBuffer_Scalar, MPI::DOUBLE,
                         Buffer_Receive_Coord_x, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_Coord_y, nBuffer_Scalar, MPI::DOUBLE,
                         Buffer_Receive_Coord_y, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  if (nDim == 3) MPI::COMM_WORLD.Gather(Buffer_Send_Coord_z, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Receive_Coord_z, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_Press, nBuffer_Scalar, MPI::DOUBLE,
                         Buffer_Receive_Press, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_CPress, nBuffer_Scalar, MPI::DOUBLE,
                         Buffer_Receive_CPress, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  if (config->GetKind_Solver() == EULER)
    MPI::COMM_WORLD.Gather(Buffer_Send_Mach, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Receive_Mach, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  if ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver() == RANS))
    MPI::COMM_WORLD.Gather(Buffer_Send_SkinFriction, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Receive_SkinFriction, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  
  MPI::COMM_WORLD.Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                         Buffer_Receive_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG, MASTER_NODE);
  
  /*--- The master node writes the surface files ---*/
  if (rank == MASTER_NODE) {
    
    /*--- Write file name with extension if unsteady ---*/
    char buffer[50];
    string filename = config->GetSurfFlowCoeff_FileName();
    
    /*--- Remove the domain number from the surface csv filename ---*/
    if (nProcessor > 1) filename.erase (filename.end()-2, filename.end());
    
    /*--- Write file name with extension if unsteady ---*/
    strcpy (cstr, filename.c_str());
    if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
      if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.csv", int(val_iZone));
      if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.csv", int(val_iZone));
      if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.csv", int(val_iZone));
      if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.csv", int(val_iZone));
      if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.csv", int(val_iZone));
      
    } else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
      if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.csv", int(iExtIter));
      if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.csv",  int(iExtIter));
      if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.csv",   int(iExtIter));
      if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv",    int(iExtIter));
      if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
    }
    else
      sprintf (buffer, ".csv");
    
    strcat (cstr, buffer);
    SurfFlow_file.precision(15);
    SurfFlow_file.open(cstr, ios::out);
    
    SurfFlow_file << "\"Global_Index\", \"x_coord\", \"y_coord\", ";
    if (nDim == 3) SurfFlow_file << "\"z_coord\", ";
    SurfFlow_file << "\"Pressure\", \"Pressure_Coefficient\", ";
    
    switch (solver) {
      case EULER : SurfFlow_file <<  "\"Mach_Number\"" << endl; break;
      case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"Skin_Friction_Coefficient\"" << endl; break;
    }
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          Global_Index = Buffer_Receive_GlobalIndex[position];
          xCoord = Buffer_Receive_Coord_x[position];
          yCoord = Buffer_Receive_Coord_y[position];
          if (nDim == 3) zCoord = Buffer_Receive_Coord_z[position];
          Pressure = Buffer_Receive_Press[position];
          PressCoeff = Buffer_Receive_CPress[position];
          SurfFlow_file << scientific << Global_Index << ", " << xCoord << ", " << yCoord << ", ";
          if (nDim == 3) SurfFlow_file << scientific << zCoord << ", ";
          SurfFlow_file << scientific << Pressure << ", " << PressCoeff << ", ";
          switch (solver) {
            case EULER :
              Mach = Buffer_Receive_Mach[position];
              SurfFlow_file << scientific << Mach << endl;
              break;
            case NAVIER_STOKES: case RANS:
              SkinFrictionCoeff = Buffer_Receive_SkinFriction[position];
              SurfFlow_file << scientific << SkinFrictionCoeff << endl;
              break;
          }
        }
      }
    }

  }
  
  if (rank == MASTER_NODE) {
    delete [] Buffer_Receive_Coord_x;
    delete [] Buffer_Receive_Coord_y;
    if (nDim == 3) delete [] Buffer_Receive_Coord_z;
    delete [] Buffer_Receive_Press;
    delete [] Buffer_Receive_CPress;
    delete [] Buffer_Receive_Mach;
    delete [] Buffer_Receive_SkinFriction;
    delete [] Buffer_Receive_GlobalIndex;
  }
  
  delete [] Buffer_Send_Coord_x;
  delete [] Buffer_Send_Coord_y;
  delete [] Buffer_Send_Coord_z;
  delete [] Buffer_Send_Press;
  delete [] Buffer_Send_CPress;
  delete [] Buffer_Send_Mach;
  delete [] Buffer_Send_SkinFriction;
  delete [] Buffer_Send_GlobalIndex;
  
  SurfFlow_file.close();
  
#endif
}

void COutput::SetSurfaceCSV_Adjoint(CConfig *config, CGeometry *geometry, CSolver *AdjSolver, CSolver *FlowSolution, unsigned long iExtIter, unsigned short val_iZone) {
  
#ifdef NO_MPI
  
  unsigned long iPoint, iVertex;
  double *Solution, xCoord, yCoord, zCoord, *IntBoundary_Jump;
  unsigned short iMarker;
  char cstr[200], buffer[50];
  ofstream SurfAdj_file;
  
  /*--- Write file name with extension if unsteady ---*/
  strcpy (cstr, config->GetSurfAdjCoeff_FileName().c_str());
  
  if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.csv", int(val_iZone));
    if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.csv", int(val_iZone));
    if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.csv", int(val_iZone));
    if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.csv", int(val_iZone));
    if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.csv", int(val_iZone));
    
  } else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.csv", int(iExtIter));
    if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.csv",  int(iExtIter));
    if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.csv",   int(iExtIter));
    if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv",    int(iExtIter));
    if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
  }
  else
    sprintf (buffer, ".csv");
  
  strcat(cstr, buffer);
  SurfAdj_file.precision(15);
  SurfAdj_file.open(cstr, ios::out);
  
  if (geometry->GetnDim() == 2) {
    SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"" << endl;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Solution = AdjSolver->node[iPoint]->GetSolution();
          IntBoundary_Jump = AdjSolver->node[iPoint]->GetIntBoundary_Jump();
          xCoord = geometry->node[iPoint]->GetCoord(0);
          yCoord = geometry->node[iPoint]->GetCoord(1);
          SurfAdj_file << scientific << iPoint << ", " << AdjSolver->GetCSensitivity(iMarker,iVertex) << ", " << Solution[0] << ", "
          << Solution[1] << ", " << Solution[2] << ", " << Solution[3] <<", " << xCoord <<", "<< yCoord << endl;
        }
    }
  }
  
  if (geometry->GetnDim() == 3) {
    SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"" << endl;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Solution = AdjSolver->node[iPoint]->GetSolution();
          xCoord = geometry->node[iPoint]->GetCoord(0);
          yCoord = geometry->node[iPoint]->GetCoord(1);
          zCoord = geometry->node[iPoint]->GetCoord(2);
          
          SurfAdj_file << scientific << iPoint << ", " << AdjSolver->GetCSensitivity(iMarker,iVertex) << ", " << Solution[0] << ", "
          << Solution[1] << ", " << Solution[2] << ", " << Solution[3] << ", " << Solution[4] << ", "<< xCoord <<", "<< yCoord <<", "<< zCoord << endl;
        }
    }
  }
  
  SurfAdj_file.close();
  
#else
  
  int rank = MPI::COMM_WORLD.Get_rank(), iProcessor, nProcessor = MPI::COMM_WORLD.Get_size();
  unsigned short nDim = geometry->GetnDim(), iMarker;
  double *Solution, *Normal, *d, *Coord;
  unsigned long Buffer_Send_nVertex[1], iVertex, iPoint, nVertex_Surface = 0, nLocalVertex_Surface = 0,
  MaxLocalVertex_Surface = 0, nBuffer_Scalar;
  unsigned long *Buffer_Receive_nVertex = NULL;
  ofstream SurfAdj_file;
  
  /*--- Write the surface .csv file ---*/
  nLocalVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface ++;
      }
  
  if (rank == MASTER_NODE)
    Buffer_Receive_nVertex = new unsigned long [nProcessor];
  
  Buffer_Send_nVertex[0] = nLocalVertex_Surface;
  
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI::UNSIGNED_LONG, MPI::MAX);
  MPI::COMM_WORLD.Gather(&Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
  
  double *Buffer_Send_Coord_x = new double[MaxLocalVertex_Surface];
  double *Buffer_Send_Coord_y= new double[MaxLocalVertex_Surface];
  double *Buffer_Send_Coord_z= new double[MaxLocalVertex_Surface];
  unsigned long *Buffer_Send_GlobalPoint= new unsigned long[MaxLocalVertex_Surface];
  double *Buffer_Send_Sensitivity= new double[MaxLocalVertex_Surface];
  double *Buffer_Send_PsiRho= new double[MaxLocalVertex_Surface];
  double *Buffer_Send_Phi_x= new double[MaxLocalVertex_Surface];
  double *Buffer_Send_Phi_y= new double[MaxLocalVertex_Surface];
  double *Buffer_Send_Phi_z= new double[MaxLocalVertex_Surface];
  double *Buffer_Send_PsiE= new double[MaxLocalVertex_Surface];
  
  nVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          Solution = AdjSolver->node[iPoint]->GetSolution();
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->node[iPoint]->GetCoord();
          d = AdjSolver->node[iPoint]->GetForceProj_Vector();
          Buffer_Send_GlobalPoint[nVertex_Surface] = geometry->node[iPoint]->GetGlobalIndex();
          Buffer_Send_Coord_x[nVertex_Surface] = Coord[0];
          Buffer_Send_Coord_y[nVertex_Surface] = Coord[1];
          Buffer_Send_Sensitivity[nVertex_Surface] =  AdjSolver->GetCSensitivity(iMarker,iVertex);
          Buffer_Send_PsiRho[nVertex_Surface] = Solution[0];
          Buffer_Send_Phi_x[nVertex_Surface] = Solution[1];
          Buffer_Send_Phi_y[nVertex_Surface] = Solution[2];
          if (nDim == 2) Buffer_Send_PsiE[nVertex_Surface] = Solution[3];
          if (nDim == 3) {
            Buffer_Send_Coord_z[nVertex_Surface] = Coord[2];
            Buffer_Send_Phi_z[nVertex_Surface] = Solution[3];
            Buffer_Send_PsiE[nVertex_Surface] = Solution[4];
          }
          nVertex_Surface++;
        }
      }
  
  double *Buffer_Receive_Coord_x = NULL, *Buffer_Receive_Coord_y = NULL, *Buffer_Receive_Coord_z = NULL, *Buffer_Receive_Sensitivity = NULL,
  *Buffer_Receive_PsiRho = NULL, *Buffer_Receive_Phi_x = NULL, *Buffer_Receive_Phi_y = NULL, *Buffer_Receive_Phi_z = NULL,
  *Buffer_Receive_PsiE = NULL;
  unsigned long *Buffer_Receive_GlobalPoint = NULL;
  
  if (rank == MASTER_NODE) {
    Buffer_Receive_Coord_x = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Coord_y = new double [nProcessor*MaxLocalVertex_Surface];
    if (nDim == 3) Buffer_Receive_Coord_z = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Sensitivity = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_PsiRho = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Phi_x = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Phi_y = new double [nProcessor*MaxLocalVertex_Surface];
    if (nDim == 3) Buffer_Receive_Phi_z = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_PsiE = new double [nProcessor*MaxLocalVertex_Surface];
  }
  
  nBuffer_Scalar = MaxLocalVertex_Surface;
  
  /*--- Send the information to the Master node ---*/
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Gather(Buffer_Send_Coord_x, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_x, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_Coord_y, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_y, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  if (nDim == 3) MPI::COMM_WORLD.Gather(Buffer_Send_Coord_z, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_z, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_GlobalPoint, nBuffer_Scalar, MPI::UNSIGNED_LONG, Buffer_Receive_GlobalPoint, nBuffer_Scalar, MPI::UNSIGNED_LONG, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_Sensitivity, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Sensitivity, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_PsiRho, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_PsiRho, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_Phi_x, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Phi_x, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_Phi_y, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Phi_y, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  if (nDim == 3) MPI::COMM_WORLD.Gather(Buffer_Send_Phi_z, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Phi_z, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_PsiE, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_PsiE, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
  
  /*--- The master node is the one who writes the surface files ---*/
  if (rank == MASTER_NODE) {
    unsigned long iVertex, GlobalPoint, position;
    char cstr[200], buffer[50];
    ofstream SurfAdj_file;
    string filename = config->GetSurfAdjCoeff_FileName();
    
    /*--- Remove the domain number from the surface csv filename ---*/
    if (nProcessor > 1) filename.erase (filename.end()-2, filename.end());
    
    /*--- Write file name with extension if unsteady ---*/
    strcpy (cstr, filename.c_str());
    
    if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
      if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.csv", int(val_iZone));
      if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.csv", int(val_iZone));
      if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.csv", int(val_iZone));
      if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.csv", int(val_iZone));
      if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.csv", int(val_iZone));
      
    } else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
      if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.csv", int(iExtIter));
      if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.csv", int(iExtIter));
      if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.csv", int(iExtIter));
      if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv", int(iExtIter));
      if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
    }
    else
      sprintf (buffer, ".csv");
    
    strcat (cstr, buffer);
    SurfAdj_file.open(cstr, ios::out);
    SurfAdj_file.precision(15);
    
    /*--- Write the 2D surface flow coefficient file ---*/
    if (geometry->GetnDim() == 2) {
      
      SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"" << endl;
      
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
        for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
          
          position = iProcessor*MaxLocalVertex_Surface+iVertex;
          GlobalPoint = Buffer_Receive_GlobalPoint[position];
          
          SurfAdj_file << scientific << GlobalPoint <<
          ", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
          ", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] <<
          ", " << Buffer_Receive_PsiE[position] << ", " << Buffer_Receive_Coord_x[position] <<
          ", "<< Buffer_Receive_Coord_y[position]  << endl;
        }
    }
    
    /*--- Write the 3D surface flow coefficient file ---*/
    if (geometry->GetnDim() == 3) {
      
      SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"" << endl;
      
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
        for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
          position = iProcessor*MaxLocalVertex_Surface+iVertex;
          GlobalPoint = Buffer_Receive_GlobalPoint[position];
          
          SurfAdj_file << scientific << GlobalPoint <<
          ", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
          ", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] << ", " << Buffer_Receive_Phi_z[position] <<
          ", " << Buffer_Receive_PsiE[position] <<", "<< Buffer_Receive_Coord_x[position] <<
          ", "<< Buffer_Receive_Coord_y[position] <<", "<< Buffer_Receive_Coord_z[position] << endl;
        }
    }
    
  }
  
  if (rank == MASTER_NODE) {
    delete [] Buffer_Receive_nVertex;
    delete [] Buffer_Receive_Coord_x;
    delete [] Buffer_Receive_Coord_y;
    if (nDim == 3) delete [] Buffer_Receive_Coord_z;
    delete [] Buffer_Receive_Sensitivity;
    delete [] Buffer_Receive_PsiRho;
    delete [] Buffer_Receive_Phi_x;
    delete [] Buffer_Receive_Phi_y;
    if (nDim == 3) delete [] Buffer_Receive_Phi_z;
    delete [] Buffer_Receive_PsiE;
    delete [] Buffer_Receive_GlobalPoint;
  }
  
  delete [] Buffer_Send_Coord_x;
  delete [] Buffer_Send_Coord_y;
  delete [] Buffer_Send_Coord_z;
  delete [] Buffer_Send_GlobalPoint;
  delete [] Buffer_Send_Sensitivity;
  delete [] Buffer_Send_PsiRho;
  delete [] Buffer_Send_Phi_x;
  delete [] Buffer_Send_Phi_y;
  delete [] Buffer_Send_Phi_z;
  delete [] Buffer_Send_PsiE;
  
  SurfAdj_file.close();
  
#endif
}

void COutput::SetSurfaceCSV_Linearized(CConfig *config, CGeometry *geometry, CSolver *LinSolution, string val_filename, unsigned long iExtIter) { }

void COutput::MergeConnectivity(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifndef NO_MPI
  rank = MPI::COMM_WORLD.Get_rank();
  size = MPI::COMM_WORLD.Get_size();
#endif
  
  /*--- Merge connectivity for each type of element (excluding halos). Note
   that we only need to merge the connectivity once, as it does not change
   during computation. Check whether the base file has been written. ---*/
  
  if (!wrote_base_file) {
    
    /*--- Merge volumetric grid. ---*/
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Tria != 0))
      cout <<"Merging volumetric triangle grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, TRIANGLE    );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Quad != 0))
      cout <<"Merging volumetric rectangle grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, RECTANGLE   );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Tetr != 0))
      cout <<"Merging volumetric tetrahedron grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, TETRAHEDRON );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Hexa != 0))
      cout <<"Merging volumetric hexahedron grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, HEXAHEDRON  );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Wedg != 0))
      cout <<"Merging volumetric wedge grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, WEDGE       );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Pyra != 0))
      cout <<"Merging volumetric pyramid grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, PYRAMID     );
    
    /*--- Merge surface grid. ---*/
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Line != 0))
      cout <<"Merging surface line grid connectivity." << endl;
    MergeSurfaceConnectivity(config, geometry, LINE);
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_BoundTria != 0))
      cout <<"Merging surface triangle grid connectivity." << endl;
    MergeSurfaceConnectivity(config, geometry, TRIANGLE);
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_BoundQuad != 0))
      cout <<"Merging surface rectangle grid connectivity." << endl;
    MergeSurfaceConnectivity(config, geometry, RECTANGLE);
    
    /*--- Update total number of volume elements after merge. ---*/
    
    nGlobal_Elem = nGlobal_Tria + nGlobal_Quad + nGlobal_Tetr +
    nGlobal_Hexa + nGlobal_Pyra + nGlobal_Wedg;
    
    /*--- Update total number of surface elements after merge. ---*/
    
    nSurf_Elem = nGlobal_Line + nGlobal_BoundTria + nGlobal_BoundQuad;
    
    /*--- Write the connectivity to the base binary output file, then
     clear the memory immediately for the rest of the computation. ---*/
    
    unsigned short FileFormat = config->GetOutput_FileFormat();
    if (rank == MASTER_NODE && FileFormat == CGNS_SOL) {
      SetCGNS_Connectivity(config, geometry, val_iZone);
      DeallocateConnectivity(config, geometry, false);
    }
    
  }
}

void COutput::MergeCoordinates(CConfig *config, CGeometry *geometry) {
  
  /*--- Local variables needed on all processors ---*/
  
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, jPoint;
  
#ifdef NO_MPI
  
  /*--- In serial, the single process has access to all geometry, so simply
   load the coordinates into the data structure. ---*/
  
  /*--- Total number of points in the mesh (excluding halos). ---*/
  nGlobal_Poin = geometry->GetnPointDomain();
  nGlobal_Doma = geometry->GetnPointDomain();
  
  /*--- Allocate the coordinates data structure. ---*/
  
  Coords = new double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Coords[iDim] = new double[nGlobal_Poin];
  }
  
  /*--- Loop over the mesh to collect the coords of the local points. ---*/
  
  jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Retrieve the current coordinates at this node. ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        Coords[iDim][jPoint] = geometry->node[iPoint]->GetCoord(iDim);
      }
      
      /*--- Increment a counter since we may be skipping over
       some halo nodes during this loop. ---*/
      jPoint++;
    }
  }
  
#else
  
  /*--- MPI preprocessing ---*/
  int iProcessor;
  int nProcessor = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
  
  bool Wrt_Halo = config->GetWrt_Halo();
  
  /*--- Local variables needed for merging the geometry with MPI. ---*/
  
  unsigned long Buffer_Send_nPoin[1], *Buffer_Recv_nPoin = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  
  if (rank == MASTER_NODE) Buffer_Recv_nPoin = new unsigned long[nProcessor];
  
  /*--- Sum total number of nodes that belong to the domain ---*/
  
  Buffer_Send_nPoin[0] = geometry->GetnPointDomain();
  MPI::COMM_WORLD.Gather(&Buffer_Send_nPoin, 1, MPI::UNSIGNED_LONG,
                         Buffer_Recv_nPoin, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
  if (rank == MASTER_NODE) {
    nGlobal_Doma = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Doma += Buffer_Recv_nPoin[iProcessor];
    }
  }
  
  /*--- Each processor sends its local number of nodes to the master. ---*/
  
  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else
    nLocalPoint = geometry->GetnPointDomain();
  Buffer_Send_nPoin[0] = nLocalPoint;
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allreduce(&nLocalPoint, &MaxLocalPoint,
                            1, MPI::UNSIGNED_LONG, MPI::MAX);
  MPI::COMM_WORLD.Gather(&Buffer_Send_nPoin, 1, MPI::UNSIGNED_LONG,
                         Buffer_Recv_nPoin, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
  nBuffer_Scalar = MaxLocalPoint;
  
  /*--- Send and Recv buffers. ---*/
  
  double *Buffer_Send_X = new double[MaxLocalPoint];
  double *Buffer_Recv_X = NULL;
  
  double *Buffer_Send_Y = new double[MaxLocalPoint];
  double *Buffer_Recv_Y = NULL;
  
  double *Buffer_Send_Z, *Buffer_Recv_Z = NULL;
  if (nDim == 3) Buffer_Send_Z = new double[MaxLocalPoint];
  
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Prepare the receive buffers in the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    
    Buffer_Recv_X = new double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Y = new double[nProcessor*MaxLocalPoint];
    if (nDim == 3) Buffer_Recv_Z = new double[nProcessor*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];
    
    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoin[iProcessor];
    }
    Coords = new double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Coords[iDim] = new double[nGlobal_Poin];
    }
  }
  
  /*--- Main communication routine. Loop over each coordinate and perform
   the MPI comm. Temporary 1-D buffers are used to send the coordinates at
   all nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/
  
  /*--- Loop over this partition to collect the coords of the local points.
   Note that we are NOT including the halo nodes here. ---*/
  double *Coords_Local; jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check for halos and write only if requested ---*/
    if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
      
      /*--- Retrieve local coordinates at this node. ---*/
      Coords_Local = geometry->node[iPoint]->GetCoord();
      
      /*--- Load local coords into the temporary send buffer. ---*/
      Buffer_Send_X[jPoint] = Coords_Local[0];
      Buffer_Send_Y[jPoint] = Coords_Local[1];
      if (nDim == 3) Buffer_Send_Z[jPoint] = Coords_Local[2];
      
      /*--- Store the global index for this local node. ---*/
      Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
      
      /*--- Increment jPoint as the counter. We need this because iPoint
       may include halo nodes that we skip over during this loop. ---*/
      jPoint++;
    }
  }
  
  /*--- Gather the coordinate data on the master node using MPI. ---*/
  
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Gather(Buffer_Send_X, nBuffer_Scalar, MPI::DOUBLE,
                         Buffer_Recv_X, nBuffer_Scalar, MPI::DOUBLE,
                         MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_Y, nBuffer_Scalar, MPI::DOUBLE,
                         Buffer_Recv_Y, nBuffer_Scalar, MPI::DOUBLE,
                         MASTER_NODE);
  if (nDim == 3) {
    MPI::COMM_WORLD.Gather(Buffer_Send_Z, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Z, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
  }
  MPI::COMM_WORLD.Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                         Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                         MASTER_NODE);
  
  /*--- The master node unpacks and sorts this variable by global index ---*/
  
  if (rank == MASTER_NODE) {
    jPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        
        /*--- Get global index, then loop over each variable and store ---*/
        iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
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
  
  delete [] Buffer_Send_X;
  delete [] Buffer_Send_Y;
  if (nDim == 3) delete [] Buffer_Send_Z;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_X;
    delete [] Buffer_Recv_Y;
    if (nDim == 3)  delete [] Buffer_Recv_Z;
    delete [] Buffer_Recv_GlobalIndex;
    delete [] Buffer_Recv_nPoin;
  }
  
#endif
  
}

void COutput::MergeVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
  rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- Local variables needed on all processors ---*/
  
  unsigned short NODES_PER_ELEMENT;
  
  unsigned long iPoint, iNode, jNode;
  unsigned long iElem = 0, jElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;
  
  int *Conn_Elem;
  
  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/
  
  switch (Elem_Type) {
    case TRIANGLE:
      nLocalElem = geometry->GetnElemTria();
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case RECTANGLE:
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
    case WEDGE:
      nLocalElem = geometry->GetnElemWedg();
      NODES_PER_ELEMENT = N_POINTS_WEDGE;
      break;
    case PYRAMID:
      nLocalElem = geometry->GetnElemPyra();
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      cout << "Error: Unrecognized element type \n";
      exit(0); break;
  }
  
  /*--- Merge the connectivity in serial or parallel. ---*/
  
#ifdef NO_MPI
  
  /*--- In serial, the single process has access to all connectivity,
   so simply load it into the data structure. ---*/
  
  /*--- Allocate a temporary array for the connectivity ---*/
  Conn_Elem = new int[nLocalElem*NODES_PER_ELEMENT];
  
  /*--- Load all elements of the current type into the buffer
   to be sent to the master node. ---*/
  jNode = 0; jElem = 0; nElem_Total = 0; bool isHalo;
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if(geometry->elem[iElem]->GetVTK_Type() == Elem_Type) {
      
      /*--- Check if this is a halo node. ---*/
      isHalo = false;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        iPoint = geometry->elem[iElem]->GetNode(iNode);
        if (!geometry->node[iPoint]->GetDomain())
          isHalo = true;
      }
      
      /*--- Loop over all nodes in this element and load the
       connectivity into the temporary array. Do not merge any
       halo cells (periodic BC). Note that we are adding one to
       the index value because CGNS/Tecplot use 1-based indexing. ---*/
      
      if (!isHalo) {
        nElem_Total++;
        for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
          Conn_Elem[jNode] = (int)geometry->elem[iElem]->GetNode(iNode) + 1;
          
          /*--- Increment jNode as the counter. ---*/
          jNode++;
        }
      }
    }
  }
  
#else
  
  /*--- MPI preprocessing ---*/
  
  int iProcessor, jProcessor;
  int nProcessor = MPI::COMM_WORLD.Get_size();
  
  /*--- Local variables needed for merging the geometry with MPI. ---*/
  
  unsigned long iVertex, iMarker;
  
  int SendRecv, RecvFrom;
  
  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0, pElem = 0;
  unsigned long MaxLocalElem = 0;
  
  bool Wrt_Halo = config->GetWrt_Halo();
  bool *Write_Elem;
  
  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/
  
  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[nProcessor];
  
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allreduce(&nLocalElem, &MaxLocalElem,
                            1, MPI::UNSIGNED_LONG, MPI::MAX);
  MPI::COMM_WORLD.Gather(&Buffer_Send_nElem, 1, MPI::UNSIGNED_LONG,
                         Buffer_Recv_nElem, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
  
  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;
  
  /*--- Send and Recv buffers ---*/
  
  int *Buffer_Send_Elem = new int[nBuffer_Scalar];
  int *Buffer_Recv_Elem = NULL;
  
  int *Buffer_Send_Halo = new int[MaxLocalElem];
  int *Buffer_Recv_Halo = NULL;
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = false;
  
  /*--- Prepare the receive buffers on the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new int[nProcessor*nBuffer_Scalar];
    Buffer_Recv_Halo = new int[nProcessor*MaxLocalElem];
    Conn_Elem = new int[nProcessor*MaxLocalElem*NODES_PER_ELEMENT];
  }
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the lower rank processor. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      if (SendRecv < 0 && RecvFrom < rank) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Local_Halo[iPoint] = true;
        }
      }
    }
  }
  
  /*--- Loop over all elements in this partition and load the
   elements of the current type into the buffer to be sent to
   the master node. ---*/
  
  jNode = 0; jElem = 0;
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if(geometry->elem[iElem]->GetVTK_Type() == Elem_Type) {
      
      /*--- Loop over all nodes in this element and load the
       connectivity into the send buffer. ---*/
      
      Buffer_Send_Halo[jElem] = false;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        
        /*--- Store the global index values directly. ---*/
        
        iPoint = geometry->elem[iElem]->GetNode(iNode);
        Buffer_Send_Elem[jNode] = (int)geometry->node[iPoint]->GetGlobalIndex();
        
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
  
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI::INT,
                         Buffer_Recv_Elem, nBuffer_Scalar, MPI::INT,
                         MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_Halo, MaxLocalElem, MPI::INT,
                         Buffer_Recv_Halo, MaxLocalElem, MPI::INT,
                         MASTER_NODE);
  
  /*--- The master node unpacks and sorts the connectivity. ---*/
  
  if (rank == MASTER_NODE) {
    
    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/
    
    Write_Elem = new bool[nProcessor*MaxLocalElem];
    for (iElem = 0; iElem < nProcessor*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }
    
    /*--- Remove the rind layer from the solution only if requested ---*/
    
    if (!Wrt_Halo) {
      
      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/
      
      kElem = 0;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
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
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
        
        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {
          
          /*--- Increment total count for this element type ---*/
          nElem_Total++;
          
          /*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/
          
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
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
  delete [] Local_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }
  
#endif
  
  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/
  
  if (rank == MASTER_NODE) {
    switch (Elem_Type) {
      case TRIANGLE:
        nGlobal_Tria = nElem_Total;
        if (nGlobal_Tria > 0) Conn_Tria = Conn_Elem;
        break;
      case RECTANGLE:
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
      case WEDGE:
        nGlobal_Wedg = nElem_Total;
        if (nGlobal_Wedg > 0) Conn_Wedg = Conn_Elem;
        break;
      case PYRAMID:
        nGlobal_Pyra = nElem_Total;
        if (nGlobal_Pyra > 0) Conn_Pyra = Conn_Elem;
        break;
      default:
        cout << "Error: Unrecognized element type \n";
        exit(0); break;
    }
  }
  
}

void COutput::MergeSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
  rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- Local variables needed on all processors ---*/
  
  unsigned short NODES_PER_ELEMENT;
  
  unsigned short iMarker;
  unsigned long iPoint, iNode, jNode;
  unsigned long iElem = 0, jElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;
  
  int *Conn_Elem;
  
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
    case RECTANGLE:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      cout << "Error: Unrecognized element type \n";
      exit(0); break;
  }
  
  /*--- Merge the connectivity in serial or parallel. ---*/
  
#ifdef NO_MPI
  
  /*--- In serial, the single process has access to all connectivity,
   so simply load it into the data structure. ---*/
  
  /*--- Allocate a temporary array for the connectivity ---*/
  Conn_Elem = new int[nLocalElem*NODES_PER_ELEMENT];
  
  /*--- Load all elements of the current type into the buffer
   to be sent to the master node. ---*/
  jNode = 0; jElem = 0; nElem_Total = 0; bool isHalo;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
        
        if (geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
          
          /*--- Check if this is a halo node. ---*/
          isHalo = false;
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
            if (!geometry->node[iPoint]->GetDomain())
              isHalo = true;
          }
          
          /*--- Loop over all nodes in this element and load the
           connectivity into the temporary array. Do not merge any
           halo cells (periodic BC). Note that we are adding one to
           the index value because CGNS/Tecplot use 1-based indexing. ---*/
          if (!isHalo) {
            nElem_Total++;
            for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
              Conn_Elem[jNode] = (int)geometry->bound[iMarker][iElem]->GetNode(iNode) + 1;
              
              /*--- Increment jNode as the counter. ---*/
              jNode++;
            }
          }
        }
      }
  
#else
  
  /*--- MPI preprocessing ---*/
  
  int iProcessor, jProcessor;
  int nProcessor = MPI::COMM_WORLD.Get_size();
  
  /*--- Local variables needed for merging the geometry with MPI. ---*/
  
  unsigned long iVertex;
  
  int SendRecv, RecvFrom;
  
  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0, pElem = 0;
  unsigned long MaxLocalElem = 0;
  
  bool Wrt_Halo = config->GetWrt_Halo();
  bool *Write_Elem;
  
  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/
  
  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[nProcessor];
  
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allreduce(&nLocalElem, &MaxLocalElem,
                            1, MPI::UNSIGNED_LONG, MPI::MAX);
  MPI::COMM_WORLD.Gather(&Buffer_Send_nElem, 1, MPI::UNSIGNED_LONG,
                         Buffer_Recv_nElem, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
  
  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;
  
  /*--- Send and Recv buffers ---*/
  
  int *Buffer_Send_Elem = new int[nBuffer_Scalar];
  int *Buffer_Recv_Elem = NULL;
  
  int *Buffer_Send_Halo = new int[MaxLocalElem];
  int *Buffer_Recv_Halo = NULL;
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = false;
  
  /*--- Prepare the receive buffers on the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new int[nProcessor*nBuffer_Scalar];
    Buffer_Recv_Halo = new int[nProcessor*MaxLocalElem];
    Conn_Elem = new int[nProcessor*MaxLocalElem*NODES_PER_ELEMENT];
  }
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the lower rank processor. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      if (SendRecv < 0 && RecvFrom < rank) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Local_Halo[iPoint] = true;
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
      for(iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
        
        if(geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
          
          /*--- Loop over all nodes in this element and load the
           connectivity into the send buffer. ---*/
          
          Buffer_Send_Halo[jElem] = false;
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            
            /*--- Store the global index values directly. ---*/
            
            iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
            Buffer_Send_Elem[jNode] = (int)geometry->node[iPoint]->GetGlobalIndex();
            
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
  
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI::INT,
                         Buffer_Recv_Elem, nBuffer_Scalar, MPI::INT,
                         MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_Halo, MaxLocalElem, MPI::INT,
                         Buffer_Recv_Halo, MaxLocalElem, MPI::INT,
                         MASTER_NODE);
  
  /*--- The master node unpacks and sorts the connectivity. ---*/
  
  if (rank == MASTER_NODE) {
    
    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/
    
    Write_Elem = new bool[nProcessor*MaxLocalElem];
    for (iElem = 0; iElem < nProcessor*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }
    
    /*--- Remove the rind layer from the solution only if requested ---*/
    
    if (!Wrt_Halo) {
      
      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/
      
      kElem = 0;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
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
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
        
        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {
          
          /*--- Increment total count for this element type ---*/
          nElem_Total++;
          
          /*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/
          
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
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
  delete [] Local_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }
  
#endif
  
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
      case RECTANGLE:
        nGlobal_BoundQuad = nElem_Total;
        if (nGlobal_BoundQuad > 0) Conn_BoundQuad = Conn_Elem;
        break;
      default:
        cout << "Error: Unrecognized element type \n";
        exit(0); break;
    }
  }
  
}

void COutput::MergeSolution(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  /*--- Local variables needed on all processors ---*/
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, jVar, iSpecies, FirstIndex = NONE, SecondIndex = NONE, ThirdIndex = NONE;
  unsigned short nVar_First = 0, nVar_Second = 0, nVar_Third = 0, iVar_Eddy = 0, iVar_Sharp = 0;
  unsigned short iVar_GridVel = 0, iVar_PressMach = 0, iVar_Density = 0, iVar_TempLam = 0,
  iVar_Tempv = 0,iVar_MagF = 0, iVar_EF =0, iVar_Temp = 0, iVar_Lam =0, iVar_Mach = 0, iVar_Press = 0,
  iVar_ViscCoeffs = 0, iVar_Sens = 0, iVar_FEA = 0, iVar_Extra = 0;
  
  unsigned long iPoint = 0, jPoint = 0, iVertex = 0, iMarker = 0;
  double Gas_Constant, Mach2Vel, Mach_Motion, RefDensity, RefPressure, factor;

  double *Aux_Frict, *Aux_Heat, *Aux_yPlus, *Aux_Sens;
  
  bool grid_movement  = (config->GetGrid_Movement());
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface    = (config->GetKind_Regime() == FREESURFACE);
  bool transition     = (config->GetKind_Trans_Model() == LM);
  bool flow           = (( config->GetKind_Solver() == EULER             ) ||
                         ( config->GetKind_Solver() == NAVIER_STOKES     ) ||
                         ( config->GetKind_Solver() == RANS              ) ||
                         ( config->GetKind_Solver() == ADJ_EULER         ) ||
                         ( config->GetKind_Solver() == ADJ_NAVIER_STOKES ) ||
                         ( config->GetKind_Solver() == ADJ_RANS          )   );
  
  unsigned short iDim;
  bool nDim               = geometry->GetnDim();
  double RefAreaCoeff     = config->GetRefAreaCoeff();
  double Gamma            = config->GetGamma();
  double RefVel2;
  
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
    factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
  }
  
  /*--- Prepare send buffers for the conservative variables. Need to
   find the total number of conservative variables and also the
   index for their particular solution container. ---*/
  switch (Kind_Solver) {
    case EULER : case NAVIER_STOKES:
      FirstIndex = FLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case RANS :
      FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL;
      if (transition) ThirdIndex=TRANS_SOL;
      else ThirdIndex = NONE;
      break;
    case TNE2_EULER :
      FirstIndex = TNE2_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case POISSON_EQUATION:
      FirstIndex = POISSON_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case WAVE_EQUATION:
      FirstIndex = WAVE_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case HEAT_EQUATION:
      FirstIndex = HEAT_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case LINEAR_ELASTICITY:
      FirstIndex = FEA_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case ADJ_EULER : case ADJ_NAVIER_STOKES :
      FirstIndex = ADJFLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case ADJ_TNE2_EULER : case ADJ_TNE2_NAVIER_STOKES :
      FirstIndex = ADJTNE2_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case ADJ_RANS :
      FirstIndex = ADJFLOW_SOL;
      if (config->GetFrozen_Visc()) SecondIndex = NONE;
      else SecondIndex = ADJTURB_SOL;
      ThirdIndex = NONE;
      break;
    case LIN_EULER : case LIN_NAVIER_STOKES : ThirdIndex = NONE;
      FirstIndex = LINFLOW_SOL; SecondIndex = NONE;
      break;
    default: SecondIndex = NONE; ThirdIndex = NONE;
      break;
  }
  nVar_First = solver[FirstIndex]->GetnVar();
  if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
  if (ThirdIndex != NONE) nVar_Third = solver[ThirdIndex]->GetnVar();
  nVar_Consv = nVar_First + nVar_Second + nVar_Third;
  
  if (config->GetWrt_Residuals()) nVar_Total = 2*nVar_Consv;
  else nVar_Total = nVar_Consv;
  
  /*--- Add the grid velocity to the restart file for the unsteady adjoint ---*/
  if (grid_movement) {
    iVar_GridVel = nVar_Total;
    if (geometry->GetnDim() == 2) nVar_Total += 2;
    else if (geometry->GetnDim() == 3) nVar_Total += 3;
  }
  
  if ((config->GetKind_Regime() == FREESURFACE)) {
    /*--- Density ---*/
    iVar_Density = nVar_Total;
    nVar_Total += 1;
  }
  
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    /*--- Pressure, Cp, Mach ---*/
    iVar_PressMach = nVar_Total;
    nVar_Total += 3;
  }
  
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    /*--- Temperature, Laminar Viscosity ---*/
    iVar_TempLam = nVar_Total;
    nVar_Total += 2;
    /*--- Skin Friction, Heat Flux, & yPlus ---*/
    iVar_ViscCoeffs = nVar_Total;
    nVar_Total += 3;
  }
  
  if (Kind_Solver == RANS) {
    /*--- Eddy Viscosity ---*/
    iVar_Eddy = nVar_Total;
    nVar_Total += 1;
  }
  
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    /*--- Sharp edges ---*/
    iVar_Sharp = nVar_Total;
    nVar_Total += 1;
  }
  
  if ((Kind_Solver == TNE2_EULER)         ||
      (Kind_Solver == TNE2_NAVIER_STOKES)   ) {
    /*--- Mach ---*/
    iVar_Mach = nVar_Total;
    nVar_Total++;
    /*--- Pressure ---*/
    iVar_Press = nVar_Total;
    nVar_Total++;
    /*--- Temperature ---*/
    iVar_Temp = nVar_Total;
    nVar_Total++;
    /*--- Vib-El. Temperature ---*/
    iVar_Tempv = nVar_Total;
    nVar_Total++;
  }
  
  if (Kind_Solver == TNE2_NAVIER_STOKES) {
    /*--- Diffusivity, viscosity, & thermal conductivity ---*/
    iVar_TempLam = nVar_Total;
    nVar_Total += config->GetnSpecies()+3;
  }
  
  if (Kind_Solver == POISSON_EQUATION) {
    iVar_EF = geometry->GetnDim();
    nVar_Total += geometry->GetnDim();
  }
  
  if (( Kind_Solver == ADJ_EULER              ) ||
      ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
      ( Kind_Solver == ADJ_RANS               ) ||
      ( Kind_Solver == ADJ_TNE2_EULER         ) ||
      ( Kind_Solver == ADJ_TNE2_NAVIER_STOKES )   ) {
    /*--- Surface sensitivity coefficient, and solution sensor ---*/
    iVar_Sens   = nVar_Total;
    nVar_Total += 2;
  }
  
  if (Kind_Solver == LINEAR_ELASTICITY) {
    /*--- Surface sensitivity coefficient, and solution sensor ---*/
    iVar_FEA   = nVar_Total;
    nVar_Total += 2;
  }
  
  if (config->GetExtraOutput()) {
    if (Kind_Solver == RANS) {
      iVar_Extra  = nVar_Total;
      nVar_Extra  = solver[TURB_SOL]->GetnOutputVariables();
      nVar_Total += nVar_Extra;
    }
    if ((Kind_Solver == TNE2_EULER)         ||
        (Kind_Solver == TNE2_NAVIER_STOKES)) {
      iVar_Extra  = nVar_Total;
      nVar_Extra  = solver[TNE2_SOL]->GetnVar();
      nVar_Total += nVar_Extra;
    }
  }
  
  /*--- Merge the solution either in serial or parallel. ---*/
  
#ifdef NO_MPI
  
  /*--- In serial, the single process has access to all solution data,
   so it is simple to retrieve and store inside Solution_Data. ---*/
  nGlobal_Poin = geometry->GetnPointDomain();
  Data = new double*[nVar_Total];
  for (iVar = 0; iVar < nVar_Total; iVar++) {
    Data[iVar] = new double[nGlobal_Poin];
  }
  
  /*--- In case there is grid movement ---*/
  double *Grid_Vel;
  
  /*--- First, loop through the mesh in order to find and store the
   value of some coefficients at any surface nodes. They
   will be placed in an auxiliary vector and then communicated like
   all other volumetric variables. ---*/
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    
    Aux_Frict = new double [geometry->GetnPointDomain()];
    Aux_Heat  = new double [geometry->GetnPointDomain()];
    Aux_yPlus = new double [geometry->GetnPointDomain()];
    
    for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
      Aux_Frict[iPoint] = 0.0;
      Aux_Heat[iPoint]  = 0.0;
      Aux_yPlus[iPoint] = 0.0;
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Frict[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker,iVertex);
          Aux_Heat[iPoint]  = solver[FLOW_SOL]->GetHeatTransferCoeff(iMarker,iVertex);
          Aux_yPlus[iPoint] = solver[FLOW_SOL]->GetYPlus(iMarker,iVertex);
        }
      }
  }
  
  if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) ||
      (Kind_Solver == ADJ_RANS)  || (Kind_Solver == ADJ_TNE2_EULER)    ||
      (Kind_Solver == ADJ_TNE2_NAVIER_STOKES)) {
    
    Aux_Sens = new double [geometry->GetnPointDomain()];
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Aux_Sens[iPoint] = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Sens[iPoint] = solver[ADJFLOW_SOL]->GetCSensitivity(iMarker,iVertex);
        }
      }
    
  }
  
  /*--- Loop over all points in the mesh, but only write data
   for nodes in the domain (ignore periodic halo nodes). ---*/
  jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check for halo nodes & only write if requested ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Solution (first, second and third system of equations) ---*/
      jVar = 0;
      for (iVar = 0; iVar < nVar_First; iVar++) {
        Data[jVar][jPoint] = solver[FirstIndex]->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
      
      for (iVar = 0; iVar < nVar_Second; iVar++) {
        Data[jVar][jPoint] = solver[SecondIndex]->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
      
      for (iVar = 0; iVar < nVar_Third; iVar++) {
        Data[jVar][jPoint] = solver[ThirdIndex]->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
      
      /*--- Residual (first, second and third system of equations) ---*/
      if (config->GetWrt_Residuals()) {
        for (iVar = 0; iVar < nVar_First; iVar++) {
          Data[jVar][jPoint] = solver[FirstIndex]->LinSysRes.GetBlock(iPoint, iVar);
          jVar++;
        }
        
        for (iVar = 0; iVar < nVar_Second; iVar++) {
          Data[jVar][jPoint] = solver[SecondIndex]->LinSysRes.GetBlock(iPoint, iVar);
          jVar++;
        }
        
        for (iVar = 0; iVar < nVar_Third; iVar++) {
          Data[jVar][jPoint] = solver[ThirdIndex]->LinSysRes.GetBlock(iPoint, iVar);
          jVar++;
        }
      }
      
      /*--- For unsteady problems with grid movement, write the mesh velocities ---*/
      if (grid_movement) {
        Grid_Vel = geometry->node[iPoint]->GetGridVel();
        for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
          Data[jVar][jPoint] = Grid_Vel[iDim];
          jVar++;
        }
      }
      
      /*--- Any extra data for output files ---*/
      switch (Kind_Solver) {
        case EULER:

          /*--- Load buffers with the pressure, Cp, and mach variables. ---*/
          if (compressible) {
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(COMPRESSIBLE); jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure(COMPRESSIBLE) - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
            Data[jVar][jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance(); jVar++;
          }
          if (incompressible || freesurface) {
            if (freesurface) {
              Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc(); jVar++;
            }
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(INCOMPRESSIBLE); jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure(INCOMPRESSIBLE) - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref())); jVar++;
            Data[jVar][jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance(); jVar++;
            
          }
          break;
          /*--- Write pressure, Cp, mach, temperature, laminar viscosity, skin friction, heat transfer, yplus ---*/
        case NAVIER_STOKES:
          if (compressible) {
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(COMPRESSIBLE); jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure(COMPRESSIBLE) - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/
            solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
            Data[jVar][jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(); jVar++;
            Data[jVar][jPoint] = Aux_Frict[iPoint]; jVar++;
            Data[jVar][jPoint] = Aux_Heat[iPoint];  jVar++;
            Data[jVar][jPoint] = Aux_yPlus[iPoint]; jVar++;
          }
          if (incompressible || freesurface) {
            if (freesurface) {
              Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc(); jVar++;
            }
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(INCOMPRESSIBLE); jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure(INCOMPRESSIBLE) - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref())); jVar++;
            Data[jVar][jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance(); jVar++;
            Data[jVar][jPoint] = 0.0; jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(); jVar++;
            Data[jVar][jPoint] = Aux_Frict[iPoint]; jVar++;
            Data[jVar][jPoint] = Aux_Heat[iPoint];  jVar++;
            Data[jVar][jPoint] = Aux_yPlus[iPoint]; jVar++;
          }
          break;
          /*--- Write pressure, Cp, mach, temperature, laminar viscosity, skin friction, heat transfer, yplus, eddy viscosity ---*/
        case RANS:
          if (compressible) {
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(COMPRESSIBLE); jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure(COMPRESSIBLE) - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/
            solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(); jVar++;
            Data[jVar][jPoint] = Aux_Frict[iPoint]; jVar++;
            Data[jVar][jPoint] = Aux_Heat[iPoint];  jVar++;
            Data[jVar][jPoint] = Aux_yPlus[iPoint]; jVar++;
          }
          if (incompressible || freesurface) {
            if (freesurface) {
              Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc(); jVar++;
            }
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(INCOMPRESSIBLE); jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure(INCOMPRESSIBLE) - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref())); jVar++;
            Data[jVar][jPoint] = 0.0; jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(); jVar++;
            Data[jVar][jPoint] = Aux_Frict[iPoint]; jVar++;
            Data[jVar][jPoint] = Aux_Heat[iPoint];  jVar++;
            Data[jVar][jPoint] = Aux_yPlus[iPoint]; jVar++;
          }
          Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosity(); jVar++;
          Data[jVar][jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance(); jVar++;
          break;
          /*--- Write poisson field. ---*/
        case POISSON_EQUATION:
          for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
            Data[jVar][jPoint] = -1.0*solver[POISSON_SOL]->node[iPoint]->GetGradient(0,iDim);
            jVar++;
          }
          break;
          
        case TNE2_EULER:
          /*--- Write Mach number ---*/
          Data[jVar][jPoint] = sqrt(solver[TNE2_SOL]->node[iPoint]->GetVelocity2())
          / solver[TNE2_SOL]->node[iPoint]->GetSoundSpeed();
          jVar++;
          /*--- Write Pressure ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetPressure();
          jVar++;
          /*--- Write Temperature ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetTemperature();
          jVar++;
          /*--- Write Vib.-El. Temperature ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetTemperature_ve();
          jVar++;
          break;
          
        case TNE2_NAVIER_STOKES:
          /*--- Write Mach number ---*/
          Data[jVar][jPoint] = sqrt(solver[TNE2_SOL]->node[iPoint]->GetVelocity2())
          / solver[TNE2_SOL]->node[iPoint]->GetSoundSpeed();
          jVar++;
          /*--- Write Pressure ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetPressure();
          jVar++;
          /*--- Write Temperature ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetTemperature();
          jVar++;
          /*--- Write Vib.-El. Temperature ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetTemperature_ve();
          jVar++;
          /*--- Write species diffusion coefficients ---*/
          for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
            Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetDiffusionCoeff()[iSpecies];
            jVar++;
          }
          /*--- Write viscosity ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetLaminarViscosity();
          jVar++;
          /*--- Write thermal conductivity ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetThermalConductivity();
          jVar++;
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetThermalConductivity_ve();
          break;
          
        case ADJ_EULER:      case ADJ_NAVIER_STOKES:     case ADJ_RANS:
        case ADJ_TNE2_EULER: case ADJ_TNE2_NAVIER_STOKES:
          
          Data[jVar][jPoint] = Aux_Sens[iPoint]; jVar++;
          if (config->GetKind_ConvNumScheme() == SPACE_CENTERED)
          { Data[jVar][jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetSensor(); jVar++; }
          if (config->GetKind_ConvNumScheme() == SPACE_UPWIND)
          { Data[jVar][jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetLimiter(0); jVar++; }
          break;
          
        case LINEAR_ELASTICITY:
          Data[jVar][jPoint] = solver[FEA_SOL]->node[iPoint]->GetVonMises_Stress();
          jVar++;
          Data[jVar][jPoint] = solver[FEA_SOL]->node[iPoint]->GetFlow_Pressure();
          jVar++;
          break;
          
      }
    }
    
    if (config->GetExtraOutput()) {
      if (Kind_Solver == RANS) {
        for (unsigned short iVar = 0; iVar < nVar_Extra; iVar++) {
          Data[jVar][jPoint] =  solver[TURB_SOL]->OutputVariables[iPoint*nVar_Extra+iVar];
          jVar++;
        }
      }
      if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES)) {
        for (unsigned short iVar = 0; iVar < nVar_Extra; iVar++) {
          Data[jVar][jPoint] =  solver[TNE2_SOL]->OutputVariables[iPoint*nVar_Extra+iVar];
          jVar++;
        }
      }
    }
    
    /*--- Increment jPoint as the counter. We need this because iPoint
     may include halo nodes that we skip over during this loop. ---*/
    jPoint++;
    
  }
  
#else
  
  /*--- MPI preprocessing ---*/
  
  int rank = MPI::COMM_WORLD.Get_rank();
  int nProcessor = MPI::COMM_WORLD.Get_size();
  int iProcessor;
  
  /*--- Local variables needed for merging with MPI ---*/
  unsigned short CurrentIndex;
  
  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  
  bool Wrt_Halo = config->GetWrt_Halo();
  
  /*--- Each processor sends its local number of nodes to the master. ---*/
  
  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else
    nLocalPoint = geometry->GetnPointDomain();
  Buffer_Send_nPoint[0] = nLocalPoint;
  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI::UNSIGNED_LONG, MPI::MAX);
  MPI::COMM_WORLD.Gather(&Buffer_Send_nPoint, 1, MPI::UNSIGNED_LONG,
                         Buffer_Recv_nPoint, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
  nBuffer_Scalar = MaxLocalPoint;
  
  /*--- Send and Recv buffers. ---*/
  
  double *Buffer_Send_Var = new double[MaxLocalPoint];
  double *Buffer_Recv_Var = NULL;
  
  double *Buffer_Send_Res = new double[MaxLocalPoint];
  double *Buffer_Recv_Res = NULL;
  
  double *Buffer_Send_Vol = new double[MaxLocalPoint];
  double *Buffer_Recv_Vol = NULL;
  
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Auxiliary vectors for surface coefficients ---*/
  
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    Aux_Frict = new double[geometry->GetnPoint()];
    Aux_Heat  = new double[geometry->GetnPoint()];
    Aux_yPlus = new double[geometry->GetnPoint()];
  }
  
  if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) || (Kind_Solver == ADJ_TNE2_EULER) || (Kind_Solver == ADJ_TNE2_NAVIER_STOKES)) {
    Aux_Sens = new double[geometry->GetnPoint()];
    
  }
  
  /*--- Prepare the receive buffers in the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    
    Buffer_Recv_Var = new double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Res = new double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Vol = new double[nProcessor*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];
    
    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
    }
    Data = new double*[nVar_Total];
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      Data[iVar] = new double[nGlobal_Poin];
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
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Get this variable into the temporary send buffer. ---*/
        Buffer_Send_Var[jPoint] = solver[CurrentIndex]->node[iPoint]->GetSolution(jVar);
        if (config->GetWrt_Residuals()) {
          Buffer_Send_Res[jPoint] = solver[CurrentIndex]->LinSysRes.GetBlock(iPoint, jVar);
        }
        
        /*--- Only send/recv the volumes & global indices during the first loop ---*/
        if (iVar == 0) {
          Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
        }
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    if (config->GetWrt_Residuals()) {
      MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                             Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                             MASTER_NODE);
    }
    if (iVar == 0) {
      MPI::COMM_WORLD.Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                             Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                             MASTER_NODE);
    }
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          if (config->GetWrt_Residuals()) {
            Data[iVar+nVar_Consv][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          }
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  /*--- Additional communication routine for the grid velocity. Note that
   we are reusing the same temporary buffers from above for efficiency.
   Also, in the future more routines like this could be used to write
   an arbitrary number of additional variables to the file. ---*/
  
  if (grid_movement) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0; double *Grid_Vel;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the three grid velocity components. ---*/
        Grid_Vel = geometry->node[iPoint]->GetGridVel();
        Buffer_Send_Var[jPoint] = Grid_Vel[0];
        Buffer_Send_Res[jPoint] = Grid_Vel[1];
        if (geometry->GetnDim() == 3) Buffer_Send_Vol[jPoint] = Grid_Vel[2];
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    if (geometry->GetnDim() == 3) {
      MPI::COMM_WORLD.Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI::DOUBLE,
                             Buffer_Recv_Vol, nBuffer_Scalar, MPI::DOUBLE,
                             MASTER_NODE);
    }
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_GridVel;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
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
  
  /*--- Communicate the Density in Free-surface problems ---*/
  if (config->GetKind_Regime() == FREESURFACE) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the pressure and mach variables. ---*/
        Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_Density;
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
  
  /*--- Communicate Pressure, Cp, and Mach ---*/
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    
    /*--- First, loop through the mesh in order to find and store the
     value of the coefficient of pressure at any surface nodes. They
     will be placed in an auxiliary vector and then communicated like
     all other volumetric variables. ---*/
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the pressure, Cp, and mach variables. ---*/
        if (compressible) {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(COMPRESSIBLE);
          Buffer_Send_Res[jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure(COMPRESSIBLE) - RefPressure)*factor*RefAreaCoeff;
          Buffer_Send_Vol[jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/
          solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed();
        }
        if (incompressible || freesurface) {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(INCOMPRESSIBLE);
          Buffer_Send_Res[jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure(INCOMPRESSIBLE) - RefPressure)*factor*RefAreaCoeff;
          Buffer_Send_Vol[jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/
          sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref()));
        }
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    MPI::COMM_WORLD.Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Vol, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_PressMach;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
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
  
  /*--- Communicate Temperature & Laminar Viscosity ---*/
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
        if (compressible) {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature();
          Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
        }
        if (incompressible || freesurface) {
          Buffer_Send_Var[jPoint] = 0.0;
          Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc();
        }
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_TempLam;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
          Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  /*--- Communicate skin friction, heat transfer, y+ ---*/
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    
    /*--- First, loop through the mesh in order to find and store the
     value of the viscous coefficients at any surface nodes. They
     will be placed in an auxiliary vector and then communicated like
     all other volumetric variables. ---*/
    
    for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      Aux_Frict[iPoint] = 0.0;
      Aux_Heat[iPoint]  = 0.0;
      Aux_yPlus[iPoint] = 0.0;
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Frict[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker,iVertex);
          Aux_Heat[iPoint]  = solver[FLOW_SOL]->GetHeatTransferCoeff(iMarker,iVertex);
          Aux_yPlus[iPoint] = solver[FLOW_SOL]->GetYPlus(iMarker,iVertex);
        }
      }
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/
        if (compressible) {
          Buffer_Send_Var[jPoint] = Aux_Frict[iPoint];
          Buffer_Send_Res[jPoint] = Aux_Heat[iPoint];
          Buffer_Send_Vol[jPoint] = Aux_yPlus[iPoint];
        }
        if (incompressible || freesurface) {
          Buffer_Send_Var[jPoint] = Aux_Frict[iPoint];
          Buffer_Send_Res[jPoint] = Aux_Heat[iPoint];
          Buffer_Send_Vol[jPoint] = Aux_yPlus[iPoint];
        }
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    MPI::COMM_WORLD.Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Vol, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_ViscCoeffs;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar+0][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  /*--- Communicate the Eddy Viscosity ---*/
  if (Kind_Solver == RANS) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the pressure and mach variables. ---*/
        if (compressible) {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosity();
        }
        if (incompressible || freesurface) {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosity();
        }
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_Eddy;
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
  
  /*--- Communicate the Sharp Edges ---*/
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the pressure and mach variables. ---*/
        Buffer_Send_Var[jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_Sharp;
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
  
  /*--- Communicate additional variables for the two-temperature solvers ---*/
  if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES)) {
    
    /*--- Mach number ---*/
    // Loop over this partition to collect the current variable
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the Mach number variables. ---*/
        Buffer_Send_Var[jPoint] =
        sqrt(solver[TNE2_SOL]->node[iPoint]->GetVelocity2())
        /solver[TNE2_SOL]->node[iPoint]->GetSoundSpeed();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      iVar = iVar_Mach;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
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
    
    /*--- Pressure ---*/
    // Loop over this partition to collect the current variable
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the Mach number variables. ---*/
        Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetPressure();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      iVar = iVar_Press;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
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
    
    /*--- Temperature ---*/
    // Loop over this partition to collect the current variable
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the Mach number variables. ---*/
        Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetTemperature();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      iVar = iVar_Temp;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index             = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
    
    /*--- Vib-el Temperature ---*/
    // Loop over this partition to collect the current variable
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the Mach number variables. ---*/
        Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetTemperature_ve();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      iVar = iVar_Tempv;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index             = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  if (Kind_Solver == TNE2_NAVIER_STOKES) {
    /*--- Species diffusion coefficients ---*/
    // Loop over this partition to collect the current variable
    for (unsigned short iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
          
          /*--- Load buffers with the Mach number variables. ---*/
          Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetDiffusionCoeff()[iSpecies];
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      MPI::COMM_WORLD.Barrier();
      MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                             Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                             MASTER_NODE);
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      if (rank == MASTER_NODE) {
        jPoint = 0;
        iVar = iVar_TempLam+iSpecies;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
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
    
    /*--- Laminar viscosity ---*/
    // Loop over this partition to collect the current variable
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the Mach number variables. ---*/
        Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetLaminarViscosity();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      iVar = iVar_TempLam+config->GetnSpecies();
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
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
    
    /*--- Thermal conductivity ---*/
    // Loop over this partition to collect the current variable
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the Mach number variables. ---*/
        Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetThermalConductivity();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      iVar = iVar_TempLam+config->GetnSpecies()+1;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
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
    
    /*--- Vib-el Thermal conductivity ---*/
    // Loop over this partition to collect the current variable
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the Mach number variables. ---*/
        Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetThermalConductivity_ve();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      iVar = iVar_TempLam+config->GetnSpecies()+2;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
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
  
  /*--- Communicate the surface sensitivity ---*/
  
  if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) || (Kind_Solver == ADJ_TNE2_EULER) || (Kind_Solver == ADJ_TNE2_NAVIER_STOKES)) {
    
    /*--- First, loop through the mesh in order to find and store the
     value of the surface sensitivity at any surface nodes. They
     will be placed in an auxiliary vector and then communicated like
     all other volumetric variables. ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Aux_Sens[iPoint] = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Sens[iPoint] = solver[ADJFLOW_SOL]->GetCSensitivity(iMarker,iVertex);
        }
      }
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/
        
        Buffer_Send_Var[jPoint] = Aux_Sens[iPoint];
        if (config->GetKind_ConvNumScheme() == SPACE_CENTERED)
          Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetSensor(iPoint);
        if (config->GetKind_ConvNumScheme() == SPACE_UPWIND)
          Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetLimiter(0);
        
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_Sens;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar+0][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  /*--- Communicate the Linear elasticity ---*/
  if ( Kind_Solver == LINEAR_ELASTICITY ) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
        Buffer_Send_Var[jPoint] = solver[FEA_SOL]->node[iPoint]->GetVonMises_Stress();
        Buffer_Send_Res[jPoint] = solver[FEA_SOL]->node[iPoint]->GetFlow_Pressure();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_FEA;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
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
        
        if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
          
          /*--- Get this variable into the temporary send buffer. ---*/
          if (Kind_Solver == RANS) {
            Buffer_Send_Var[jPoint] = solver[TURB_SOL]->OutputVariables[iPoint*nVar_Extra+jVar];
          }
          if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES)) {
            Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->OutputVariables[iPoint*nVar_Extra+jVar];
          }
          jPoint++;
          
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      MPI::COMM_WORLD.Barrier();
      MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                             Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                             MASTER_NODE);
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_Extra;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
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
  
  
  /*--- Immediately release the temporary buffers. ---*/
  
  delete [] Buffer_Send_Var;
  delete [] Buffer_Send_Res;
  delete [] Buffer_Send_Vol;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_Var;
    delete [] Buffer_Recv_Res;
    delete [] Buffer_Recv_Vol;
    delete [] Buffer_Recv_GlobalIndex;
  }
  
#endif
  
  /*--- Release memory needed for surface coefficients ---*/
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    delete [] Aux_Frict; delete [] Aux_Heat; delete [] Aux_yPlus;
  }
  if (( Kind_Solver == ADJ_EULER              ) ||
      ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
      ( Kind_Solver == ADJ_RANS               ) ||
      ( Kind_Solver == ADJ_TNE2_EULER         ) ||
      ( Kind_Solver == ADJ_TNE2_NAVIER_STOKES )   ) {
    delete [] Aux_Sens;
  }
  
}

void COutput::MergeBaselineSolution(CConfig *config, CGeometry *geometry, CSolver *solver, unsigned short val_iZone) {
  
  /*--- Local variables needed on all processors ---*/
  unsigned short iVar;
  unsigned long iPoint = 0, jPoint = 0;
  
  nVar_Total = config->fields.size() - 1;
  
  /*--- Merge the solution either in serial or parallel. ---*/
#ifdef NO_MPI
  
  /*--- In serial, the single process has access to all solution data,
   so it is simple to retrieve and store inside Solution_Data. ---*/
  nGlobal_Poin = geometry->GetnPointDomain();
  Data = new double*[nVar_Total];
  for (iVar = 0; iVar < nVar_Total; iVar++) {
    Data[iVar] = new double[nGlobal_Poin];
  }
  
  /*--- Loop over all points in the mesh, but only write data
   for nodes in the domain (ignore periodic halo nodes). ---*/
  jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    if (geometry->node[iPoint]->GetDomain()) {
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
  
  int rank = MPI::COMM_WORLD.Get_rank();
  int nProcessor = MPI::COMM_WORLD.Get_size();
  int iProcessor;
  
  /*--- Local variables needed for merging with MPI ---*/
  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  
  bool Wrt_Halo = config->GetWrt_Halo();
  
  /*--- Each processor sends its local number of nodes to the master. ---*/
  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else
    nLocalPoint = geometry->GetnPointDomain();
  Buffer_Send_nPoint[0] = nLocalPoint;
  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI::UNSIGNED_LONG, MPI::MAX);
  MPI::COMM_WORLD.Gather(&Buffer_Send_nPoint, 1, MPI::UNSIGNED_LONG,
                         Buffer_Recv_nPoint, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
  nBuffer_Scalar = MaxLocalPoint;
  
  /*--- Send and Recv buffers. ---*/
  
  double *Buffer_Send_Var = new double[MaxLocalPoint];
  double *Buffer_Recv_Var = NULL;
  
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Prepare the receive buffers in the master node only. ---*/
  if (rank == MASTER_NODE) {
    
    Buffer_Recv_Var = new double[nProcessor*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];
    
    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
    }
    Data = new double*[nVar_Total];
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      Data[iVar] = new double[nGlobal_Poin];
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
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
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
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    if (iVar == 0) {
      MPI::COMM_WORLD.Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                             Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                             MASTER_NODE);
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
  
}

void COutput::SetRestart(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
  
  /*--- Local variables ---*/
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iExtIter = config->GetExtIter();
  bool grid_movement = config->GetGrid_Movement();
  ofstream restart_file;
  string filename;
  
  /*--- Retrieve filename from config ---*/
  if (config->GetAdjoint()) {
    filename = config->GetRestart_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else {
    filename = config->GetRestart_FlowFileName();
  }
  
  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    filename = config->GetUnsteady_FileName(filename, int(val_iZone));
  } else if (config->GetWrt_Unsteady()) {
    filename = config->GetUnsteady_FileName(filename, int(iExtIter));
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
    restart_file << "\t\"Conservative_" << iVar+1<<"\"";
  }
  if (config->GetWrt_Residuals()) {
    for (iVar = 0; iVar < nVar_Consv; iVar++) {
      restart_file << "\t\"Residual_" << iVar+1<<"\"";
    }
  }
  
  /*--- Mesh velocities for dynamic mesh cases ---*/
  if (grid_movement) {
    if (nDim == 2) {
      restart_file << "\t\"Grid_Velx\"\t\"Grid_Vely\"";
    } else {
      restart_file << "\t\"Grid_Velx\"\t\"Grid_Vely\"\t\"Grid_Velz\"";
    }
  }
  
  /*--- Solver specific output variables ---*/
  if (config->GetKind_Regime() == FREESURFACE) {
    restart_file << "\t\"Density\"";
  }
  
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    restart_file << "\t\"Pressure\"\t\"Pressure_Coefficient\"\t\"Mach\"";
  }
  
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    restart_file << "\t\"Temperature\"\t\"Laminar_Viscosity\"\t\"Skin_Friction_Coefficient\"\t\"Heat_Transfer\"\t\"Y_Plus\"";
  }
  
  if (Kind_Solver == RANS) {
    restart_file << "\t\"Eddy_Viscosity\"";
  }
  
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    restart_file << "\t\"Sharp_Edge_Dist\"";
  }
  
  if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES)) {
    restart_file << "\t\"Mach\"\t\"Pressure\"\t\"Temperature\"\t\"Temperature_ve\"";
  }
  
  if (Kind_Solver == TNE2_NAVIER_STOKES) {
    for (unsigned short iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
      restart_file << "\t\"DiffusionCoeff_" << iSpecies << "\"";
    restart_file << "\t\"Laminar_Viscosity\"\t\"ThermConductivity\"\t\"ThermConductivity_ve\"";
  }
  
  if (Kind_Solver == POISSON_EQUATION) {
    for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
      restart_file << "\t\"poissonField_" << iDim+1 << "\"";
  }
  
  if ((Kind_Solver == ADJ_EULER              ) ||
      (Kind_Solver == ADJ_NAVIER_STOKES      ) ||
      (Kind_Solver == ADJ_RANS               ) ||
      (Kind_Solver == ADJ_TNE2_EULER         ) ||
      (Kind_Solver == ADJ_TNE2_NAVIER_STOKES )   ) {
    restart_file << "\t\"Surface_Sensitivity\"\t\"Solution_Sensor\"";
  }
  
  if (Kind_Solver == LINEAR_ELASTICITY) {
    restart_file << "\t\"Von_Mises_Stress\"\t\"Flow_Pressure\"";
  }
  
  if (config->GetExtraOutput()) {
    for (iVar = 0; iVar < nVar_Extra; iVar++) {
      restart_file << "\t\"ExtraOutput_" << iVar+1<<"\"";
    }
  }
  
  restart_file << endl;
  
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
    restart_file << endl;
  }
  
  restart_file.close();
  
}

void COutput::DeallocateCoordinates(CConfig *config, CGeometry *geometry) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
  rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- Local variables and initialization ---*/
  
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
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
  rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- The master node alone owns all data found in this routine. ---*/
  if (rank == MASTER_NODE) {
    
    /*--- Deallocate memory for connectivity data ---*/
    if (surf_sol) {
      if (nGlobal_Line > 0) delete [] Conn_Line;
      if (nGlobal_BoundTria > 0) delete [] Conn_BoundTria;
      if (nGlobal_BoundQuad > 0) delete [] Conn_BoundQuad;
    }
    else {
      if (nGlobal_Tria > 0) delete [] Conn_Tria;
      if (nGlobal_Quad > 0) delete [] Conn_Quad;
      if (nGlobal_Tetr > 0) delete [] Conn_Tetr;
      if (nGlobal_Hexa > 0) delete [] Conn_Hexa;
      if (nGlobal_Wedg > 0) delete [] Conn_Wedg;
      if (nGlobal_Pyra > 0) delete [] Conn_Pyra;
    }
    
  }
}

void COutput::DeallocateSolution(CConfig *config, CGeometry *geometry) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
  rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- The master node alone owns all data found in this routine. ---*/
  if (rank == MASTER_NODE) {
    
    /*--- Deallocate memory for solution data ---*/
    for (unsigned short iVar = 0; iVar < nVar_Total; iVar++) {
      delete [] Data[iVar];
    }
    delete [] Data;
    
  }
}

void COutput::SetHistory_Header(ofstream *ConvHist_file, CConfig *config) {
  char cstr[200], buffer[50], turb_resid[1000];
  unsigned short iMarker, iMarker_Monitoring, iSpecies;
  string Monitoring_Tag, Monitoring_coeff;
  
  bool rotating_frame = config->GetRotating_Frame();
  bool aeroelastic = config->GetAeroelastic_Simulation();
  bool equiv_area = config->GetEquivArea();
  bool turbulent = ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == ADJ_RANS));
  bool frozen_turb = config->GetFrozen_Visc();
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  
  bool isothermal = false;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) isothermal = true;
  
  /* Find the markers being monitored and create a header for them */
  /* Note that for now only the aeroelastic case will make use of this */
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
    Monitoring_coeff += ",\"plunge_" + Monitoring_Tag + "\"";
    Monitoring_coeff += ",\"pitch_"  + Monitoring_Tag + "\"";
    Monitoring_coeff += ",\"CLift_"  + Monitoring_Tag + "\"";
    Monitoring_coeff += ",\"CDrag_"  + Monitoring_Tag + "\"";
    Monitoring_coeff += ",\"CMx_"    + Monitoring_Tag + "\"";
    Monitoring_coeff += ",\"CMy_"    + Monitoring_Tag + "\"";
    Monitoring_coeff += ",\"CMz_"    + Monitoring_Tag + "\"";
  }
  
  /*--- Write file name with extension ---*/
  string filename = config->GetConv_FileName();
  strcpy (cstr, filename.data());
  
  if (config->GetWrt_Unsteady() && config->GetRestart()) {
    long iExtIter = config->GetUnst_RestartIter();
    if (int(iExtIter) < 10) sprintf (buffer, "_0000%d", int(iExtIter));
    if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d", int(iExtIter));
    if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d", int(iExtIter));
    if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d", int(iExtIter));
    if (int(iExtIter) >= 10000) sprintf (buffer, "_%d", int(iExtIter));
    strcat(cstr,buffer);
  }
  
  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == TECPLOT_BINARY))  sprintf (buffer, ".plt");
  if ((config->GetOutput_FileFormat() == CGNS_SOL) ||
      (config->GetOutput_FileFormat() == PARAVIEW))  sprintf (buffer, ".csv");
  strcat(cstr,buffer);
  
  ConvHist_file->open(cstr, ios::out);
  ConvHist_file->precision(15);
  
  /*--- Begin of the header ---*/
  
  char begin[]= "\"Iteration\"";
  
  /*--- Header for the coefficients ---*/
  
  char flow_coeff[]= ",\"CLift\",\"CDrag\",\"CSideForce\",\"CMx\",\"CMy\",\"CMz\",\"CFx\",\"CFy\",\"CFz\",\"CL/CD\"";
  char heat_coeff[]= ",\"CHeat_Load\",\"CHeat_Max\"";
  char equivalent_area_coeff[]= ",\"CEquivArea\",\"CNearFieldOF\"";
  char rotating_frame_coeff[]= ",\"CMerit\",\"CT\",\"CQ\"";
  //char aeroelastic_coeff[]= Monitoring_coeff;
  string aeroelastic_coeff = Monitoring_coeff;
  char free_surface_coeff[]= ",\"CFreeSurface\"";
  char wave_coeff[]= ",\"CWave\"";
  char fea_coeff[]= ",\"CFEA\"";
  char adj_coeff[]= ",\"Sens_Geo\",\"Sens_Mach\",\"Sens_AoA\",\"Sens_Press\",\"Sens_Temp\",\"Sens_AoS\"";
  
  /*--- Header for the residuals ---*/
  
  char flow_resid[]= ",\"Res_Flow[0]\",\"Res_Flow[1]\",\"Res_Flow[2]\",\"Res_Flow[3]\",\"Res_Flow[4]\"";
  char adj_flow_resid[]= ",\"Res_AdjFlow[0]\",\"Res_AdjFlow[1]\",\"Res_AdjFlow[2]\",\"Res_AdjFlow[3]\",\"Res_AdjFlow[4]\"";
  switch (config->GetKind_Turb_Model()) {
    case SA:	sprintf (turb_resid, ",\"Res_Turb[0]\""); break;
    case ML:	sprintf (turb_resid, ",\"Res_Turb[0]\""); break;
    case SST:	sprintf (turb_resid, ",\"Res_Turb[0]\",\"Res_Turb[1]\""); break;
  }
  char adj_turb_resid[]= ",\"Res_AdjTurb[0]\"";
  char levelset_resid[]= ",\"Res_LevelSet\"";
  char adj_levelset_resid[]= ",\"Res_AdjLevelSet\"";
  char wave_resid[]= ",\"Res_Wave[0]\",\"Res_Wave[1]\"";
  char fea_resid[]= ",\"Res_FEA\"";
  char heat_resid[]= ",\"Res_Heat\"";
  
  /*--- End of the header ---*/
  
  char end[]= ",\"Linear_Solver_Iterations\",\"Time(min)\"\n";
  
  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == TECPLOT_BINARY)) {
    ConvHist_file[0] << "TITLE = \"SU2 Simulation\"" << endl;
    ConvHist_file[0] << "VARIABLES = ";
  }
  
  /*--- Write the header, case depending ---*/
  switch (config->GetKind_Solver()) {
      
    case EULER : case NAVIER_STOKES: case RANS :
      ConvHist_file[0] << begin << flow_coeff;
      if (isothermal) ConvHist_file[0] << heat_coeff;
      if (equiv_area) ConvHist_file[0] << equivalent_area_coeff;
      if (rotating_frame) ConvHist_file[0] << rotating_frame_coeff;
      if (aeroelastic) ConvHist_file[0] << aeroelastic_coeff;
      ConvHist_file[0] << flow_resid;
      if (turbulent) ConvHist_file[0] << turb_resid;
      ConvHist_file[0] << end;
      if (freesurface) {
        ConvHist_file[0] << begin << flow_coeff << free_surface_coeff;
        ConvHist_file[0] << flow_resid << levelset_resid << end;
      }
      break;
      
    case TNE2_EULER : case TNE2_NAVIER_STOKES:
      ConvHist_file[0] << begin << flow_coeff;
      if (config->GetKind_Solver() == TNE2_NAVIER_STOKES)
        ConvHist_file[0] << heat_coeff;
      for (iSpecies = 0; iSpecies < config->GetnSpecies()+5; iSpecies++)
        ConvHist_file[0] << ",\"Residual[" << iSpecies << "]\"";
      ConvHist_file[0] << end;
      break;
      
    case ADJ_EULER      : case ADJ_NAVIER_STOKES      : case ADJ_RANS:
    case ADJ_TNE2_EULER : case ADJ_TNE2_NAVIER_STOKES :
      ConvHist_file[0] << begin << adj_coeff << adj_flow_resid;
      if ((turbulent) && (!frozen_turb)) ConvHist_file[0] << adj_turb_resid;
      ConvHist_file[0] << end;
      if (freesurface) {
        ConvHist_file[0] << begin << adj_coeff << adj_flow_resid << adj_levelset_resid << end;
      }
      break;
      
    case WAVE_EQUATION:
      ConvHist_file[0] << begin << wave_coeff;
      ConvHist_file[0] << wave_resid << end;
      break;
      
    case HEAT_EQUATION:
      ConvHist_file[0] << begin << heat_coeff;
      ConvHist_file[0] << heat_resid << end;
      break;
      
    case LINEAR_ELASTICITY:
      ConvHist_file[0] << begin << fea_coeff;
      ConvHist_file[0] << fea_resid << end;
      break;
      
  }
  
  if (config->GetOutput_FileFormat() == TECPLOT || config->GetOutput_FileFormat() == TECPLOT_BINARY) {
    ConvHist_file[0] << "ZONE T= \"Convergence history\"" << endl;
  }
  
}


void COutput::SetConvergence_History(ofstream *ConvHist_file, CGeometry ***geometry, CSolver ****solver_container, CConfig **config, CIntegration ***integration, bool DualTime_Iteration, unsigned long timeused, unsigned short val_iZone) {
  
#ifndef NO_MPI
  int rank = MPI::COMM_WORLD.Get_rank();
#else
  int rank = MASTER_NODE;
#endif
  
  /*--- Output using only the master node ---*/
  if (rank == MASTER_NODE) {
    
    unsigned long iIntIter = config[val_iZone]->GetIntIter();
    unsigned long iExtIter = config[val_iZone]->GetExtIter();
    
    /*--- WARNING: These buffers have hard-coded lengths. Note that you
     may have to adjust them to be larger if adding more entries. ---*/
    char begin[1000], direct_coeff[1000], surface_coeff[1000], adjoint_coeff[1000], flow_resid[1000], adj_flow_resid[1000],
    turb_resid[1000], trans_resid[1000], adj_turb_resid[1000], resid_aux[1000],
    levelset_resid[1000], adj_levelset_resid[1000], wave_coeff[1000], heat_coeff[1000], fea_coeff[1000], wave_resid[1000], heat_resid[1000],
    fea_resid[1000], end[1000];
    double dummy = 0.0;
    unsigned short iVar, iMarker, iMarker_Monitoring;
    unsigned short iSpecies, loc;
    
    unsigned long LinSolvIter = 0;
    double timeiter = double(timeused)/double(iExtIter+1);
    
    unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
    unsigned short nDim = geometry[val_iZone][FinestMesh]->GetnDim();
    unsigned short nSpecies = config[val_iZone]->GetnSpecies();
    
    bool compressible = (config[val_iZone]->GetKind_Regime() == COMPRESSIBLE);
    bool incompressible = (config[val_iZone]->GetKind_Regime() == INCOMPRESSIBLE);
    bool freesurface = (config[val_iZone]->GetKind_Regime() == FREESURFACE);
    
    bool rotating_frame = config[val_iZone]->GetRotating_Frame();
    bool aeroelastic = config[val_iZone]->GetAeroelastic_Simulation();
    bool equiv_area = config[val_iZone]->GetEquivArea();
    bool transition = (config[val_iZone]->GetKind_Trans_Model() == LM);
    bool isothermal = false;
    for (iMarker = 0; iMarker < config[val_iZone]->GetnMarker_All(); iMarker++)
      if (config[val_iZone]->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) isothermal = true;
    bool turbulent = ((config[val_iZone]->GetKind_Solver() == RANS) || (config[val_iZone]->GetKind_Solver() == ADJ_RANS) ||
                      (config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_RANS));
    bool adjoint = config[val_iZone]->GetAdjoint();
    bool fluid_structure = ((config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_EULER) || (config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_NAVIER_STOKES) ||
                            (config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_RANS));
    bool wave = (config[val_iZone]->GetKind_Solver() == WAVE_EQUATION);
    bool heat = (config[val_iZone]->GetKind_Solver() == HEAT_EQUATION);
    bool fea = (config[val_iZone]->GetKind_Solver() == LINEAR_ELASTICITY);
    bool TNE2 = ((config[val_iZone]->GetKind_Solver() == TNE2_EULER) || (config[val_iZone]->GetKind_Solver() == TNE2_NAVIER_STOKES) ||
                 (config[val_iZone]->GetKind_Solver() == ADJ_TNE2_EULER) || (config[val_iZone]->GetKind_Solver() == ADJ_TNE2_NAVIER_STOKES));
    bool flow = (config[val_iZone]->GetKind_Regime() == EULER) || (config[val_iZone]->GetKind_Regime() == NAVIER_STOKES) ||
    (config[val_iZone]->GetKind_Regime() == RANS) || (config[val_iZone]->GetKind_Regime() == ADJ_EULER) ||
    (config[val_iZone]->GetKind_Regime() == ADJ_NAVIER_STOKES) || (config[val_iZone]->GetKind_Regime() == ADJ_RANS);
    
    /*--- Initialize variables to store information from all domains (direct solution) ---*/
    double Total_CLift = 0.0, Total_CDrag = 0.0, Total_CSideForce = 0.0, Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0, Total_CEff = 0.0,
    Total_CEquivArea = 0.0, Total_CNearFieldOF = 0.0, Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0, Total_CMerit = 0.0,
    Total_CT = 0.0, Total_CQ = 0.0, Total_CFreeSurface = 0.0, Total_CWave = 0.0, Total_CHeat = 0.0, Total_CFEA = 0.0, PressureDrag = 0.0, ViscDrag = 0.0, MagDrag = 0.0, Total_Q = 0.0, Total_MaxQ = 0.0;
    
    /*--- Initialize variables to store information from all domains (adjoint solution) ---*/
    double Total_Sens_Geo = 0.0, Total_Sens_Mach = 0.0, Total_Sens_AoA = 0.0;
    double Total_Sens_Press = 0.0, Total_Sens_Temp = 0.0;
    
    /*--- Residual arrays ---*/
    double *residual_flow = NULL, *residual_turbulent = NULL, *residual_transition = NULL, *residual_TNE2 = NULL, *residual_levelset = NULL;
    double *residual_adjflow = NULL, *residual_adjturbulent = NULL, *residual_adjTNE2 = NULL, *residual_adjlevelset = NULL;
    double *residual_wave = NULL; double *residual_fea = NULL; double *residual_heat = NULL;
    
    /*--- Coefficients Monitored arrays ---*/
    double *aeroelastic_plunge = NULL, *aeroelastic_pitch = NULL, *Surface_CLift = NULL, *Surface_CDrag = NULL, *Surface_CMx = NULL, *Surface_CMy = NULL, *Surface_CMz = NULL;
    
    /*--- Initialize number of variables ---*/
    unsigned short nVar_Flow = 0, nVar_LevelSet = 0, nVar_Turb = 0, nVar_Trans = 0, nVar_TNE2 = 0, nVar_Wave = 0, nVar_Heat = 0, nVar_FEA = 0,     nVar_AdjFlow = 0, nVar_AdjTNE2 = 0, nVar_AdjLevelSet = 0, nVar_AdjTurb = 0;
    
    /*--- Direct problem variables ---*/
    if (compressible) nVar_Flow = nDim+2; else nVar_Flow = nDim+1;
    if (turbulent) {
      switch (config[val_iZone]->GetKind_Turb_Model()){
        case SA:	nVar_Turb = 1; break;
        case ML:	nVar_Turb = 1; break;
        case SST: nVar_Turb = 2; break;
      }
    }
    if (transition) nVar_Trans = 2;
    if (TNE2) nVar_TNE2 = config[val_iZone]->GetnSpecies()+nDim+2;
    if (wave) nVar_Wave = 2;
    if (fea) nVar_FEA = nDim;
    if (heat) nVar_Heat = 1;
    if (freesurface) nVar_LevelSet = 1;
    
    /*--- Adjoint problem variables ---*/
    if (compressible) nVar_AdjFlow = nDim+2; else nVar_AdjFlow = nDim+1;
    if (turbulent) {
      switch (config[val_iZone]->GetKind_Turb_Model()){
        case SA:	nVar_AdjTurb = 1; break;
        case ML:	nVar_AdjTurb = 1; break;
        case SST: nVar_AdjTurb = 2; break;
      }
    }
    if (TNE2) nVar_AdjTNE2 = config[val_iZone]->GetnSpecies()+nDim+2;
    if (freesurface) nVar_AdjLevelSet = 1;
    
    /*--- Allocate memory for the residual ---*/
    residual_flow       = new double[nVar_Flow];
    residual_turbulent  = new double[nVar_Turb];
    residual_transition = new double[nVar_Trans];
    residual_TNE2       = new double[nVar_TNE2];
    residual_levelset   = new double[nVar_LevelSet];
    residual_wave       = new double[nVar_Wave];
    residual_fea        = new double[nVar_FEA];
    residual_heat       = new double[nVar_Heat];
    
    residual_adjflow      = new double[nVar_AdjFlow];
    residual_adjturbulent = new double[nVar_AdjTurb];
    residual_adjTNE2      = new double[nVar_AdjTNE2];
    residual_adjlevelset  = new double[nVar_AdjLevelSet];
    
    /*--- Allocate memory for the coefficients being monitored ---*/
    aeroelastic_plunge = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    aeroelastic_pitch  = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CLift      = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CDrag      = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMx        = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMy        = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMz        = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    
    /*--- Write information from nodes ---*/
    switch (config[val_iZone]->GetKind_Solver()) {
        
      case EULER:                   case NAVIER_STOKES:                   case RANS:
      case FLUID_STRUCTURE_EULER:   case FLUID_STRUCTURE_NAVIER_STOKES:   case FLUID_STRUCTURE_RANS:
      case ADJ_EULER:               case ADJ_NAVIER_STOKES:               case ADJ_RANS:
        
        /*--- Flow solution coefficients ---*/
        Total_CLift       = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CLift();
        Total_CDrag       = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CDrag();
        Total_CSideForce  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CSideForce();
        Total_CEff        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CEff();
        Total_CMx         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMx();
        Total_CMy         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMy();
        Total_CMz         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMz();
        Total_CFx         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFx();
        Total_CFy         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFy();
        Total_CFz         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFz();
        
        if (freesurface) {
          Total_CFreeSurface = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFreeSurface();
        }
        
        if (isothermal) {
          Total_Q     = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_Q();
          Total_MaxQ  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_MaxQ();
        }
        
        if (equiv_area) {
          Total_CEquivArea    = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CEquivArea();
          Total_CNearFieldOF  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CNearFieldOF();
          
          /*--- Note that there is a redefinition of the nearfield based functionals ---*/
          Total_CEquivArea    = config[val_iZone]->GetWeightCd()*Total_CDrag + (1.0-config[val_iZone]->GetWeightCd())*Total_CEquivArea;
          Total_CNearFieldOF  = config[val_iZone]->GetWeightCd()*Total_CDrag + (1.0-config[val_iZone]->GetWeightCd())*Total_CNearFieldOF;
        }
        
        if (rotating_frame) {
          Total_CT      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CT();
          Total_CQ      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CQ();
          Total_CMerit  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMerit();
        }
        
        if (aeroelastic) {
          /*--- Look over the markers being monitored and get the desired values ---*/
          for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
            aeroelastic_plunge[iMarker_Monitoring] = config[val_iZone]->GetAeroelastic_plunge(iMarker_Monitoring);
            aeroelastic_pitch[iMarker_Monitoring]  = config[val_iZone]->GetAeroelastic_pitch(iMarker_Monitoring);
            Surface_CLift[iMarker_Monitoring]      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CLift(iMarker_Monitoring);
            Surface_CDrag[iMarker_Monitoring]      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CDrag(iMarker_Monitoring);
            Surface_CMx[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMx(iMarker_Monitoring);
            Surface_CMy[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMy(iMarker_Monitoring);
            Surface_CMz[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMz(iMarker_Monitoring);
          }
        }
        
        if (fluid_structure) {
          Total_CFEA  = solver_container[ZONE_1][FinestMesh][FEA_SOL]->GetTotal_CFEA();
        }
        
        /*--- Flow Residuals ---*/
        
        for (iVar = 0; iVar < nVar_Flow; iVar++)
          residual_flow[iVar] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_RMS(iVar);
        
        /*--- Turbulent residual ---*/
        
        if (turbulent) {
          for (iVar = 0; iVar < nVar_Turb; iVar++)
            residual_turbulent[iVar] = solver_container[val_iZone][FinestMesh][TURB_SOL]->GetRes_RMS(iVar);
        }
        
        /*--- Transition residual ---*/
        
        if (transition) {
          for (iVar = 0; iVar < nVar_Trans; iVar++)
            residual_transition[iVar] = solver_container[val_iZone][FinestMesh][TRANS_SOL]->GetRes_RMS(iVar);
        }
        
        /*--- Free Surface residual ---*/
        
        if (freesurface) {
          for (iVar = 0; iVar < nVar_LevelSet; iVar++)
            residual_levelset[iVar] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_RMS(nDim+1);
        }
        
        /*--- FEA residual ---*/
        
        if (fluid_structure) {
          for (iVar = 0; iVar < nVar_FEA; iVar++)
            residual_fea[iVar] = solver_container[ZONE_1][FinestMesh][FEA_SOL]->GetRes_RMS(iVar);
        }
        
        /*--- Iterations of the linear solver ---*/
        
        LinSolvIter = (unsigned long) solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetIterLinSolver();
        
        /*--- Adjoint solver ---*/
        
        if (adjoint) {
          
          /*--- Adjoint solution coefficients ---*/
          
          Total_Sens_Geo   = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Geo();
          Total_Sens_Mach  = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Mach();
          Total_Sens_AoA   = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_AoA();
          Total_Sens_Press = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Press();
          Total_Sens_Temp  = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Temp();
          
          /*--- Adjoint flow residuals ---*/
          
          for (iVar = 0; iVar < nVar_AdjFlow; iVar++) {
            residual_adjflow[iVar] = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(iVar);
          }
          
          /*--- Adjoint turbulent residuals ---*/
          
          if (turbulent) {
            if (!config[val_iZone]->GetFrozen_Visc()) {
              for (iVar = 0; iVar < nVar_AdjTurb; iVar++)
                residual_adjturbulent[iVar] = solver_container[val_iZone][FinestMesh][ADJTURB_SOL]->GetRes_RMS(iVar);
            }
          }
          
          /*--- Adjoint level set residuals ---*/
          
          if (freesurface) {
            for (iVar = 0; iVar < nVar_AdjLevelSet; iVar++)
              residual_adjlevelset[iVar] = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(nDim+1);
          }
          
        }
        
        break;
        
      case TNE2_EULER:     case TNE2_NAVIER_STOKES:
      case ADJ_TNE2_EULER: case ADJ_TNE2_NAVIER_STOKES:
        
        /*--- Coefficients ---*/
        
        Total_CLift       = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CLift();
        Total_CDrag       = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CDrag();
        Total_CSideForce  = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CSideForce();
        Total_CEff        = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CEff();
        Total_CMx         = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CMx();
        Total_CMy         = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CMy();
        Total_CMz         = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CMz();
        Total_CFx         = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CFx();
        Total_CFy         = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CFy();
        Total_CFz         = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CFz();
        if (config[val_iZone]->GetKind_Solver() == TNE2_NAVIER_STOKES) {
          Total_Q           = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_Q();
          Total_MaxQ        = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_MaxQ();
        }
        
        /*--- Residuals ---*/
        
        for (iVar = 0; iVar < nVar_TNE2; iVar++)
          residual_TNE2[iVar] = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetRes_RMS(iVar);
        
        /*--- Iterations of the linear solver ---*/
        LinSolvIter = (unsigned long) solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetIterLinSolver();
        
        /*--- Adjoint solver ---*/
        if (adjoint) {
          
          /*--- Adjoint solution coefficients ---*/
          
          Total_Sens_Geo   = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_Geo();
          Total_Sens_Mach  = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_Mach();
          Total_Sens_AoA   = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_AoA();
          Total_Sens_Press = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_Press();
          Total_Sens_Temp  = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_Temp();
          
          /*--- Adjoint flow residuals ---*/
          
          for (iVar = 0; iVar < nVar_AdjTNE2; iVar++) {
            residual_adjTNE2[iVar] = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetRes_RMS(iVar);
          }
        }
        
        break;
        
      case WAVE_EQUATION:
        
        /*--- Wave coefficients  ---*/
        
        Total_CWave = solver_container[val_iZone][FinestMesh][WAVE_SOL]->GetTotal_CWave();
        
        /*--- Wave Residuals ---*/
        
        for (iVar = 0; iVar < nVar_Wave; iVar++) {
          residual_wave[iVar] = solver_container[val_iZone][FinestMesh][WAVE_SOL]->GetRes_RMS(iVar);
        }
        
        break;
        
      case HEAT_EQUATION:
        
        /*--- Heat coefficients  ---*/
        
        Total_CHeat = solver_container[val_iZone][FinestMesh][HEAT_SOL]->GetTotal_CHeat();
        
        /*--- Wave Residuals ---*/
        
        for (iVar = 0; iVar < nVar_Heat; iVar++) {
          residual_heat[iVar] = solver_container[val_iZone][FinestMesh][HEAT_SOL]->GetRes_RMS(iVar);
        }
        
        break;
        
      case LINEAR_ELASTICITY:
        
        /*--- FEA coefficients ---*/
        
        Total_CFEA = solver_container[val_iZone][FinestMesh][FEA_SOL]->GetTotal_CFEA();
        
        /*--- Plasma Residuals ---*/
        
        for (iVar = 0; iVar < nVar_FEA; iVar++) {
          residual_fea[iVar] = solver_container[val_iZone][FinestMesh][FEA_SOL]->GetRes_RMS(iVar);
        }
        
        break;
        
    }
    
    /*--- Header frecuency ---*/
    
    bool Unsteady = ((config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                     (config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND));
    bool In_NoDualTime = (!DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    bool In_DualTime_0 = (DualTime_Iteration && (iIntIter % config[val_iZone]->GetWrt_Con_Freq_DualTime() == 0));
    bool In_DualTime_1 = (!DualTime_Iteration && Unsteady);
    bool In_DualTime_2 = (Unsteady && DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    bool In_DualTime_3 = (Unsteady && !DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    
    bool write_heads;
    if (Unsteady) write_heads = (iIntIter == 0);
    else write_heads = (((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*20)) == 0));
    
    if ((In_NoDualTime || In_DualTime_0 || In_DualTime_1) && (In_NoDualTime || In_DualTime_2 || In_DualTime_3)) {
      
      /*--- Prepare the history file output, note that the dual
       time output don't write to the history file ---*/
      if (!DualTime_Iteration) {
        
        /*--- Write the begining of the history file ---*/
        sprintf (begin, "%12d", int(iExtIter));
        
        /*--- Write the end of the history file ---*/
        sprintf (end, ", %12.10f, %12.10f\n", double(LinSolvIter), double(timeused)/(CLOCKS_PER_SEC*60.0));
        
        /*--- Write the solution and residual of the history file ---*/
        switch (config[val_iZone]->GetKind_Solver()) {
            
          case EULER : case NAVIER_STOKES: case RANS:
          case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES: case FLUID_STRUCTURE_RANS:
          case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
          case ADJ_TNE2_EULER: case ADJ_TNE2_NAVIER_STOKES:
            
            /*--- Direct coefficients ---*/
            sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                     Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
                     Total_CFz, Total_CEff);
            if (isothermal)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy,
                       Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_Q, Total_MaxQ);
            if (equiv_area)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy,
                       Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CEquivArea, Total_CNearFieldOF);
            if (rotating_frame)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
                       Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CMerit, Total_CT, Total_CQ);
            if (aeroelastic) {
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
                       Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff);
              for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
                //Append one by one the surface coeff to direct coeff. (Think better way do this, maybe use string)
                sprintf(surface_coeff, ", %12.10f",aeroelastic_plunge[iMarker_Monitoring]);
                strcat(direct_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",aeroelastic_pitch[iMarker_Monitoring]);
                strcat(direct_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CLift[iMarker_Monitoring]);
                strcat(direct_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CDrag[iMarker_Monitoring]);
                strcat(direct_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CMx[iMarker_Monitoring]);
                strcat(direct_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CMy[iMarker_Monitoring]);
                strcat(direct_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CMz[iMarker_Monitoring]);
                strcat(direct_coeff, surface_coeff);
              }
            }
            if (freesurface) {
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
                       Total_CFz, Total_CEff, Total_CFreeSurface);
            }
            if (fluid_structure)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz,
                       Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CFEA);
            
            /*--- Flow residual ---*/
            if (nDim == 2) {
              if (compressible) sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), dummy );
              if (incompressible || freesurface) sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), dummy, dummy );
            }
            else {
              if (compressible) sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), log10 (residual_flow[4]) );
              if (incompressible || freesurface) sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), dummy );
            }
            
            /*--- Turbulent residual ---*/
            if (turbulent){
              switch(nVar_Turb) {
                case 1: sprintf (turb_resid, ", %12.10f", log10 (residual_turbulent[0])); break;
                case 2: sprintf (turb_resid, ", %12.10f, %12.10f", log10(residual_turbulent[0]), log10(residual_turbulent[1])); break;
              }
            }
            
            /*--- Transition residual ---*/
            if (transition){
              sprintf (trans_resid, ", %12.10f, %12.10f", log10(residual_transition[0]), log10(residual_transition[1]));
            }
            
            /*--- Free surface residual ---*/
            if (freesurface) {
              sprintf (levelset_resid, ", %12.10f", log10 (residual_levelset[0]));
            }
            
            /*--- Fluid structure residual ---*/
            if (fluid_structure) {
              if (nDim == 2) sprintf (levelset_resid, ", %12.10f, %12.10f, 0.0", log10 (residual_fea[0]), log10 (residual_fea[1]));
              else sprintf (levelset_resid, ", %12.10f, %12.10f, %12.10f", log10 (residual_fea[0]), log10 (residual_fea[1]), log10 (residual_fea[2]));
            }
            
            if (adjoint) {
              
              /*--- Adjoint coefficients ---*/
              sprintf (adjoint_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, 0.0", Total_Sens_Geo, Total_Sens_Mach, Total_Sens_AoA, Total_Sens_Press, Total_Sens_Temp);
              
              /*--- Adjoint flow residuals ---*/
              if (nDim == 2) {
                if (compressible) sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0", log10 (residual_adjflow[0]),log10 (residual_adjflow[1]),log10 (residual_adjflow[2]),log10 (residual_adjflow[3]) );
                if (incompressible || freesurface) sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, 0.0, 0.0", log10 (residual_adjflow[0]),log10 (residual_adjflow[1]),log10 (residual_adjflow[2]) );
              }
              else {
                if (compressible) sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_adjflow[0]),log10 (residual_adjflow[1]),log10 (residual_adjflow[2]),log10 (residual_adjflow[3]), log10 (residual_adjflow[4]) );
                if (incompressible || freesurface) sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0", log10 (residual_adjflow[0]),log10 (residual_adjflow[1]),log10 (residual_adjflow[2]),log10 (residual_adjflow[3]) );
              }
              
              /*--- Adjoint turbulent residuals ---*/
              if (turbulent)
                if (!config[val_iZone]->GetFrozen_Visc())
                  sprintf (adj_turb_resid, ", %12.10f", log10 (residual_adjturbulent[0]));
              
              /*--- Adjoint free surface residuals ---*/
              if (freesurface) sprintf (adj_levelset_resid, ", %12.10f", log10 (residual_adjlevelset[0]));
            }
            
            break;
            
          case TNE2_EULER : case TNE2_NAVIER_STOKES:
            
            /*--- Direct coefficients ---*/
            if (config[val_iZone]->GetKind_Solver() == TNE2_NAVIER_STOKES)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                       Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
                       Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz,
                       Total_CEff, Total_Q, Total_MaxQ);
            else
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                       Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
                       Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz,
                       Total_CEff);
            
            /*--- Direct problem residual ---*/
            for (iVar = 0; iVar < nSpecies+nDim+2; iVar++) {
              sprintf (resid_aux, ", %12.10f", log10 (residual_TNE2[iVar]));
              if (iVar == 0) strcpy(flow_resid, resid_aux);
              else strcat(flow_resid, resid_aux);
            }
            
            if (adjoint) {
              
              /*--- Adjoint coefficients ---*/
              sprintf (adjoint_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, 0.0", Total_Sens_Geo, Total_Sens_Mach, Total_Sens_AoA, Total_Sens_Press, Total_Sens_Temp);
              
              /*--- Adjoint flow residuals ---*/
              for (iVar = 0; iVar < nSpecies+nDim+2; iVar++) {
                sprintf (resid_aux, ", %12.10f", log10 (residual_adjTNE2[iVar]));
                if (iVar == 0) strcpy(adj_flow_resid, resid_aux);
                else strcat(adj_flow_resid, resid_aux);
              }
            }
            
            break;
            
          case WAVE_EQUATION:
            
            sprintf (direct_coeff, ", %12.10f", Total_CWave);
            sprintf (wave_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_wave[0]), log10 (residual_wave[1]), dummy, dummy, dummy );
            
            break;
            
          case HEAT_EQUATION:
            
            sprintf (direct_coeff, ", %12.10f", Total_CHeat);
            sprintf (heat_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_heat[0]), dummy, dummy, dummy, dummy );
            
            break;
            
          case LINEAR_ELASTICITY:
            
            sprintf (direct_coeff, ", %12.10f", Total_CFEA);
            sprintf (fea_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_fea[0]), dummy, dummy, dummy, dummy );
            
            break;
            
        }
      }
      
      /*--- Write the screen header---*/
      if ((write_heads) && !(!DualTime_Iteration && Unsteady)) {
        
        if (!Unsteady) {
          switch (config[val_iZone]->GetKind_Solver()) {
            case EULER :                  case NAVIER_STOKES:
            case FLUID_STRUCTURE_EULER :  case FLUID_STRUCTURE_NAVIER_STOKES:
              cout << endl << " Min Delta Time: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<<
              ". Max Delta Time: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() << ".";
              break;
              
            case TNE2_EULER: case TNE2_NAVIER_STOKES:
            case ADJ_TNE2_EULER: case ADJ_TNE2_NAVIER_STOKES:
              cout << endl << " Min Delta Time: " << solver_container[val_iZone][MESH_0][TNE2_SOL]->GetMin_Delta_Time()<< ". Max Delta Time: " << solver_container[val_iZone][MESH_0][TNE2_SOL]->GetMax_Delta_Time() << ".";
              break;
          }
        }
        else {
          if (flow) {
            cout << endl << " Min DT: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<<
            ". Max DT: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() <<
            ". Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
          }
          else {
            cout << endl << " Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
          }
        }
        
        switch (config[val_iZone]->GetKind_Solver()) {
          case EULER :                  case NAVIER_STOKES:
          case FLUID_STRUCTURE_EULER :  case FLUID_STRUCTURE_NAVIER_STOKES:
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max(0) << "." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            
            if (!fluid_structure) {
              if (incompressible) cout << "   Res[Press]" << "     Res[Velx]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
              else if (freesurface) cout << "   Res[Press]" << "     Res[Dist]" << "   CLift(Total)" << "     CLevelSet" << endl;
              else if (rotating_frame && nDim == 3) cout << "     Res[Rho]" << "     Res[RhoE]" << " CThrust(Total)" << " CTorque(Total)" << endl;
              else if (aeroelastic) cout << "     Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << "         plunge" << "          pitch" << endl;
              else if (equiv_area) cout << "     Res[Rho]" << "   CLift(Total)" << "   CDrag(Total)" << "    CPress(N-F)" << endl;
              else cout << "     Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
            }
            else if (fluid_structure) cout << "     Res[Rho]" << "   Res[Displx]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
            
            break;
            
          case TNE2_EULER :  case TNE2_NAVIER_STOKES:
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetPoint_Max(0) << "." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            
            cout << "     Res[Rho]" << "     Res[RhoE]" << "   Res[RhoEve]" << "   CDrag(Total)";
            if (config[val_iZone]->GetKind_Solver() == TNE2_NAVIER_STOKES)
              cout << "   Max qdot" << endl;
            else cout << endl;
            break;
            
          case RANS : case FLUID_STRUCTURE_RANS:
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max(0) << "." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            if (incompressible || freesurface) cout << "   Res[Press]";
            else cout << "      Res[Rho]";
            
            switch (config[val_iZone]->GetKind_Turb_Model()){
              case SA:	cout << "       Res[nu]"; break;
              case ML:	cout << "       Res[nu]"; break;
              case SST:	cout << "     Res[kine]" << "     Res[omega]"; break;
            }
            
            if (transition) { cout << "      Res[Int]" << "       Res[Re]"; }
            if (rotating_frame && nDim == 3 ) cout << "   CThrust(Total)" << "   CTorque(Total)" << endl;
            if (aeroelastic) cout << "     plunge" << "     pitch" << endl;
            else cout << "   CLift(Total)"   << "   CDrag(Total)"   << endl;
            break;
            
          case WAVE_EQUATION :
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            cout << "      Res[Wave]" << "   CWave(Total)"<<  endl;
            break;
            
          case HEAT_EQUATION :
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            cout << "      Res[Heat]" << "   CHeat(Total)"<<  endl;
            break;
            
          case LINEAR_ELASTICITY :
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            if (nDim == 2) cout << "    Res[Displx]" << "    Res[Disply]" << "   CFEA(Total)"<<  endl;
            if (nDim == 3) cout << "    Res[Displx]" << "    Res[Disply]" << "    Res[Displz]" << "   CFEA(Total)"<<  endl;
            break;
            
          case ADJ_EULER :              case ADJ_NAVIER_STOKES :
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max(0) << "." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            if (incompressible || freesurface) cout << "   Res[Psi_Press]" << "   Res[Psi_Velx]";
            else cout << "   Res[Psi_Rho]" << "     Res[Psi_E]";
            cout << "     Sens_Geo" << "    Sens_Mach" << endl;
            
            if (freesurface) {
              cout << "   Res[Psi_Press]" << "   Res[Psi_Dist]" << "    Sens_Geo" << "   Sens_Mach" << endl;
            }
            break;
            
          case ADJ_RANS :
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max(0) << "." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            if (incompressible || freesurface) cout << "     Res[Psi_Press]";
            else cout << "     Res[Psi_Rho]";
            
            if (!config[val_iZone]->GetFrozen_Visc()) {
              cout << "      Res[Psi_nu]";
            }
            else {
              if (incompressible || freesurface) cout << "   Res[Psi_Velx]";
              else cout << "     Res[Psi_E]";
            }
            cout << "     Sens_Geo" << "    Sens_Mach" << endl;
            
            if (freesurface) {
              cout << "   Res[Psi_Press]" << "   Res[Psi_Dist]" << "    Sens_Geo" << "   Sens_Mach" << endl;
            }
            break;
            
          case ADJ_TNE2_EULER :              case ADJ_TNE2_NAVIER_STOKES :
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetPoint_Max(0) << "." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            cout << "   Res[Psi_Rho]" << "     Res[Psi_E]" << "   Res[Psi_Eve]" << "     Sens_Geo" << "    Sens_Mach" << endl;
            
            break;
            
        }
        
      }
      
      /*--- Write the solution on the screen and history file ---*/
      cout.precision(6);
      cout.setf(ios::fixed,ios::floatfield);
      
      if (!Unsteady) {
        cout.width(5); cout << iExtIter;
        cout.width(11); cout << double(timeiter)/CLOCKS_PER_SEC;
        
      } else {
        cout.width(8); cout << iIntIter;
        cout.width(8); cout << iExtIter;
      }
      
      switch (config[val_iZone]->GetKind_Solver()) {
        case EULER : case NAVIER_STOKES:
        case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES:
          
          if (!DualTime_Iteration) {
            if (compressible) ConvHist_file[0] << begin << direct_coeff << flow_resid;
            if (incompressible) ConvHist_file[0] << begin << direct_coeff << flow_resid;
            if (freesurface) ConvHist_file[0] << begin << direct_coeff << flow_resid << levelset_resid << end;
            if (fluid_structure) ConvHist_file[0] << fea_resid;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(13); cout << log10(residual_flow[0]);
          if (!fluid_structure && !equiv_area) {
            if (compressible) {
              if (nDim == 2 ) { cout.width(14); cout << log10(residual_flow[3]); }
              else { cout.width(14); cout << log10(residual_flow[4]); }
            }
            if (incompressible) { cout.width(14); cout << log10(residual_flow[1]); }
            if (freesurface) { cout.width(14); cout << log10(residual_levelset[0]); }
          }
          else if (fluid_structure) { cout.width(14); cout << log10(residual_fea[0]); }
          
          if (rotating_frame && nDim == 3 ) {
            cout.setf(ios::scientific,ios::floatfield);
            cout.width(15); cout << Total_CT;
            cout.width(15); cout << Total_CQ;
            cout.unsetf(ios_base::floatfield);
          }
          else if (equiv_area) { cout.width(15); cout << Total_CLift; cout.width(15); cout << Total_CDrag; cout.width(15);
            cout.precision(4);
            cout.setf(ios::scientific,ios::floatfield);
            cout << Total_CNearFieldOF; }
          else if (freesurface) { cout.width(15); cout << Total_CLift; cout.width(15); cout << Total_CFreeSurface; }
          else { cout.width(15); cout << min(1000.0,max(-1000.0, Total_CLift)); cout.width(15); cout << min(1000.0,max(-1000.0, Total_CDrag)); }
          if (aeroelastic) {
            cout.setf(ios::scientific,ios::floatfield);
            cout.width(15); cout << aeroelastic_plunge[0]; //Only output the first marker being monitored to the console.
            cout.width(15); cout << aeroelastic_pitch[0];
            cout.unsetf(ios_base::floatfield);
          }
          cout << endl;
          
          break;
          
        case RANS :
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << direct_coeff << flow_resid << turb_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          
          if (incompressible || freesurface) cout.width(13);
          else  cout.width(14);
          cout << log10(residual_flow[0]);
          
          switch(nVar_Turb) {
            case 1: cout.width(14); cout << log10(residual_turbulent[0]); break;
            case 2: cout.width(14); cout << log10(residual_turbulent[0]);
              cout.width(14); cout << log10(residual_turbulent[1]); break;
          }
          
          if (transition) { cout.width(14); cout << log10(residual_transition[0]); cout.width(14); cout << log10(residual_transition[1]); }
          
          if (rotating_frame  && nDim == 3 ) {
            cout.setf(ios::scientific,ios::floatfield);
            cout.width(15); cout << Total_CT; cout.width(15);
            cout << Total_CQ;
            cout.unsetf(ios_base::floatfield);
          }
          else { cout.width(15); cout << min(1000.0,max(-1000.0, Total_CLift)); cout.width(15); cout << min(1000.0,max(-1000.0, Total_CDrag)); }
          if (aeroelastic) {
            cout.setf(ios::scientific,ios::floatfield);
            cout.width(15); cout << aeroelastic_plunge[0]; //Only output the first marker being monitored to the console.
            cout.width(15); cout << aeroelastic_pitch[0];
            cout.unsetf(ios_base::floatfield);
          }
          cout << endl;
          
          if (freesurface) {
            if (!DualTime_Iteration) {
              ConvHist_file[0] << begin << direct_coeff << flow_resid << levelset_resid << end;
              ConvHist_file[0].flush();
            }
            
            cout.precision(6);
            cout.setf(ios::fixed,ios::floatfield);
            cout.width(13); cout << log10(residual_flow[0]);
            cout.width(14); cout << log10(residual_levelset[0]);
            cout.width(15); cout << Total_CLift;
            cout.width(14); cout << Total_CFreeSurface;
            
            cout << endl;
          }
          
          break;
          
        case TNE2_EULER : case TNE2_NAVIER_STOKES:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << direct_coeff << flow_resid;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(13); cout << log10(residual_TNE2[0]);
          cout.width(14); cout << log10(residual_TNE2[nSpecies+nDim]);
          cout.width(14); cout << log10(residual_TNE2[nSpecies+nDim+1]);
          cout.width(15); cout << Total_CDrag;
          if (config[val_iZone]->GetKind_Solver()==TNE2_NAVIER_STOKES) {
            cout.precision(1);
            cout.width(11); cout << Total_MaxQ;
          }
          cout << endl;
          break;
          
        case WAVE_EQUATION:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << wave_coeff << wave_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(14); cout << log10(residual_wave[0]);
          cout.width(14); cout << Total_CWave;
          cout << endl;
          break;
          
        case HEAT_EQUATION:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << heat_coeff << heat_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(14); cout << log10(residual_heat[0]);
          cout.width(14); cout << Total_CHeat;
          cout << endl;
          break;
          
        case LINEAR_ELASTICITY:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << fea_coeff << fea_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(15); cout << log10(residual_fea[0]);
          cout.width(15); cout << log10(residual_fea[1]);
          if (nDim == 3) { cout.width(15); cout << log10(residual_fea[2]); }
          cout.precision(4);
          cout.setf(ios::scientific,ios::floatfield);
          cout.width(14); cout << Total_CFEA;
          cout << endl;
          break;
          
        case ADJ_EULER :              case ADJ_NAVIER_STOKES :
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          if (compressible) {
            cout.width(15); cout << log10(residual_adjflow[0]);
            cout.width(15); cout << log10(residual_adjflow[nDim+1]);
          }
          if (incompressible || freesurface) {
            cout.width(17); cout << log10(residual_adjflow[0]);
            cout.width(16); cout << log10(residual_adjflow[1]);
          }
          cout.precision(4);
          cout.setf(ios::scientific,ios::floatfield);
          cout.width(14); cout << Total_Sens_Geo;
          cout.width(14); cout << Total_Sens_Mach;
          cout << endl;
          cout.unsetf(ios_base::floatfield);
          
          if (freesurface) {
            if (!DualTime_Iteration) {
              ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << adj_levelset_resid << end;
              ConvHist_file[0].flush();
            }
            
            cout.precision(6);
            cout.setf(ios::fixed,ios::floatfield);
            cout.width(17); cout << log10(residual_adjflow[0]);
            cout.width(16); cout << log10(residual_adjlevelset[0]);
            cout.precision(3);
            cout.setf(ios::scientific,ios::floatfield);
            cout.width(12); cout << Total_Sens_Geo;
            cout.width(12); cout << Total_Sens_Mach;
            cout.unsetf(ios_base::floatfield);
            cout << endl;
          }
          
          break;
          
        case ADJ_RANS :
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid;
            if (!config[val_iZone]->GetFrozen_Visc())
              ConvHist_file[0] << adj_turb_resid;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(17); cout << log10(residual_adjflow[0]);
          if (!config[val_iZone]->GetFrozen_Visc()) {
            cout.width(17); cout << log10(residual_adjturbulent[0]);
          }
          else {
            if (compressible) {
              if (geometry[val_iZone][FinestMesh]->GetnDim() == 2 ) { cout.width(15); cout << log10(residual_adjflow[3]); }
              else { cout.width(15); cout << log10(residual_adjflow[4]); }
            }
            if (incompressible || freesurface) {
              cout.width(15); cout << log10(residual_adjflow[1]);
            }
          }
          cout.precision(4);
          cout.setf(ios::scientific,ios::floatfield);
          cout.width(14); cout << Total_Sens_Geo;
          cout.width(14); cout << Total_Sens_Mach;
          cout << endl;
          cout.unsetf(ios_base::floatfield);
          if (freesurface) {
            if (!DualTime_Iteration) {
              ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << adj_levelset_resid;
              ConvHist_file[0] << end;
              ConvHist_file[0].flush();
            }
            
            cout.precision(6);
            cout.setf(ios::fixed,ios::floatfield);
            cout.width(17); cout << log10(residual_adjflow[0]);
            cout.width(16); cout << log10(residual_adjlevelset[0]);
            
            cout.precision(4);
            cout.setf(ios::scientific,ios::floatfield);
            cout.width(12); cout << Total_Sens_Geo;
            cout.width(12); cout << Total_Sens_Mach;
            cout << endl;
            cout.unsetf(ios_base::floatfield);
          }
          
          break;
          
        case ADJ_TNE2_EULER :              case ADJ_TNE2_NAVIER_STOKES :
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(15); cout << log10(residual_adjTNE2[0]);
          cout.width(15); cout << log10(residual_adjTNE2[nSpecies+nDim]);
          cout.width(15); cout << log10(residual_adjTNE2[nSpecies+nDim+1]);

          cout.precision(4);
          cout.setf(ios::scientific,ios::floatfield);
          cout.width(14); cout << Total_Sens_Geo;
          cout.width(14); cout << Total_Sens_Mach;
          cout << endl;
          cout.unsetf(ios_base::floatfield);
          
          break;

          
      }
      cout.unsetf(ios::fixed);
      
      delete [] residual_flow;
      delete [] residual_levelset;
      delete [] residual_TNE2;
      delete [] residual_turbulent;
      delete [] residual_transition;
      delete [] residual_wave;
      delete [] residual_fea;
      delete [] residual_heat;
      
      delete [] residual_adjflow;
      delete [] residual_adjTNE2;
      delete [] residual_adjlevelset;
      delete [] residual_adjturbulent;
      
      delete [] Surface_CLift;
      delete [] Surface_CDrag;
      delete [] Surface_CMx;
      delete [] Surface_CMy;
      delete [] Surface_CMz;
      
    }
  }
}

void COutput::SetResult_Files(CSolver ****solver_container, CGeometry ***geometry, CConfig **config,
                              unsigned long iExtIter, unsigned short val_nZone) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
  rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  unsigned short iZone;
  
  for (iZone = 0; iZone < val_nZone; iZone++) {
    
    /*--- Flags identifying the types of files to be written. ---*/
    bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
    bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();
    
#ifndef NO_MPI
    /*--- Do not merge the volume solutions if we are running in parallel.
     Force the use of SU2_SOL to merge the volume sols in this case. ---*/
    int size = MPI::COMM_WORLD.Get_size();
    if (size > SINGLE_NODE) {
      Wrt_Vol = false;
      Wrt_Srf = false;
    }
#endif
    
    bool Wrt_Csv = config[iZone]->GetWrt_Csv_Sol();
    bool Wrt_Rst = config[iZone]->GetWrt_Restart();
    
    switch (config[iZone]->GetKind_Solver()) {
        
      case EULER : case NAVIER_STOKES : case RANS :
      case FLUID_STRUCTURE_EULER : case FLUID_STRUCTURE_NAVIER_STOKES : case FLUID_STRUCTURE_RANS:
        
        if (Wrt_Csv) SetSurfaceCSV_Flow(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][FLOW_SOL], iExtIter, iZone);
        break;
        
      case TNE2_EULER : case TNE2_NAVIER_STOKES :
        
        if (Wrt_Csv) SetSurfaceCSV_Flow(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][TNE2_SOL], iExtIter, iZone);
        break;
        
      case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS :
        if (Wrt_Csv) SetSurfaceCSV_Adjoint(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][ADJFLOW_SOL], solver_container[iZone][MESH_0][FLOW_SOL], iExtIter, iZone);
        break;
        
      case LIN_EULER : case LIN_NAVIER_STOKES :
        if (Wrt_Csv) SetSurfaceCSV_Linearized(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][LINFLOW_SOL], config[iZone]->GetSurfLinCoeff_FileName(), iExtIter);
        break;
    }
    
    /*--- Get the file output format ---*/
    
    unsigned short FileFormat = config[iZone]->GetOutput_FileFormat();
    
    bool dynamic_mesh = (config[iZone]->GetUnsteady_Simulation() &&
                         config[iZone]->GetGrid_Movement());
    
    /*--- Merge the node coordinates and connectivity, if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/
    
    if (Wrt_Vol || Wrt_Srf)
      MergeConnectivity(config[iZone], geometry[iZone][MESH_0], iZone);
    
    /*--- Merge coordinates of all grid nodes (excluding ghost points).
     The grid coordinates are always merged and included first in the
     restart files. ---*/
    
    MergeCoordinates(config[iZone], geometry[iZone][MESH_0]);
    
    if (rank == MASTER_NODE) {
      
      if (FileFormat == CGNS_SOL) {
        SetCGNS_Coordinates(config[iZone], geometry[iZone][MESH_0], iZone);
        if (!wrote_base_file || dynamic_mesh)
          DeallocateCoordinates(config[iZone], geometry[iZone][MESH_0]);
      } else if (FileFormat == TECPLOT_BINARY) {
        SetTecplot_Mesh(config[iZone], geometry[iZone][MESH_0], iZone);
        SetTecplot_SurfaceMesh(config[iZone], geometry[iZone][MESH_0], iZone);
        if (!wrote_base_file)
          DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
        if (!wrote_surf_file)
          DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], wrote_surf_file);
      }
    }
    
    /*--- Merge the solution data needed for volume solutions and restarts ---*/
    
    if (Wrt_Vol || Wrt_Rst)
      MergeSolution(config[iZone], geometry[iZone][MESH_0],
                    solver_container[iZone][MESH_0], iZone);
    
    /*--- Write restart, CGNS, or Tecplot files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/
    
    if (rank == MASTER_NODE) {
      
      /*--- Write a native restart file ---*/
      if (Wrt_Rst)
        SetRestart(config[iZone], geometry[iZone][MESH_0], iZone);
      
      if (Wrt_Vol) {
        
        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            SetTecplot_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
            break;
            
          case TECPLOT_BINARY:
            
            /*--- Write a Tecplot binary solution file ---*/
            SetTecplot_Solution(config[iZone], geometry[iZone][MESH_0], iZone);
            break;
            
          case CGNS_SOL:
            
            /*--- Write a CGNS solution file ---*/
            SetCGNS_Solution(config[iZone], geometry[iZone][MESH_0], iZone);
            break;
            
          case PARAVIEW:
            
            /*--- Write a Paraview ASCII file ---*/
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
            SetTecplot_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
            break;
            
          case TECPLOT_BINARY:
            
            /*--- Write a Tecplot binary solution file ---*/
            SetTecplot_SurfaceSolution(config[iZone], geometry[iZone][MESH_0], iZone);
            break;
            
          case PARAVIEW:
            
            /*--- Write a Paraview ASCII file ---*/
            SetParaview_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
            break;
            
          default:
            break;
        }
        
      }
      
      /*--- Release memory needed for merging the solution data. ---*/
      if (((Wrt_Vol) || (Wrt_Srf)) && (FileFormat == TECPLOT ||
                                       FileFormat == TECPLOT_BINARY ||
                                       FileFormat == PARAVIEW))
        DeallocateCoordinates(config[iZone], geometry[iZone][MESH_0]);
      
      if (Wrt_Vol || Wrt_Rst)
        DeallocateSolution(config[iZone], geometry[iZone][MESH_0]);
      
    }
    
    /*--- Final broadcast (informing other procs that the base output
     file was written) & barrier to sync up after master node writes
     output files. ---*/
    
#ifndef NO_MPI
    MPI::COMM_WORLD.Bcast(&wrote_base_file, 1, MPI::INT, MASTER_NODE);
    MPI::COMM_WORLD.Barrier();
#endif
    
  }
}

void COutput::SetBaselineResult_Files(CSolver **solver, CGeometry **geometry, CConfig **config,
                                      unsigned long iExtIter, unsigned short val_nZone) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
  rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  unsigned short iZone;
  
  for (iZone = 0; iZone < val_nZone; iZone++) {
    
    /*--- Flags identifying the types of files to be written. ---*/
    
    bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
    bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();
    bool Wrt_Rst = config[iZone]->GetWrt_Restart();
    
    /*--- Get the file output format ---*/
    
    unsigned short FileFormat = config[iZone]->GetOutput_FileFormat();
    
    /*--- Merge the node coordinates and connectivity if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/
    
    if (Wrt_Vol || Wrt_Srf) {
      if (rank == MASTER_NODE) cout <<"Merging grid connectivity." << endl;
      MergeConnectivity(config[iZone], geometry[iZone], iZone);
    }
    
    /*--- Merge the solution data needed for volume solutions and restarts ---*/
    
    if (Wrt_Vol || Wrt_Rst) {
      if (rank == MASTER_NODE) cout <<"Merging solution." << endl;
      MergeBaselineSolution(config[iZone], geometry[iZone], solver[iZone], iZone);
    }
    
    /*--- Write restart, CGNS, Tecplot or Paraview files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/
    
    if (rank == MASTER_NODE) {
      
      if (Wrt_Vol) {
        
        if (rank == MASTER_NODE)
          cout <<"Writing volume solution file." << endl;
        
        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            SetTecplot_ASCII(config[iZone], geometry[iZone], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone], false);
            break;
            
          case TECPLOT_BINARY:
            
            /*--- Write a Tecplot binary solution file ---*/
            SetTecplot_Mesh(config[iZone], geometry[iZone], iZone);
            SetTecplot_Solution(config[iZone], geometry[iZone], iZone);
            break;
            
          case CGNS_SOL:
            
            /*--- Write a CGNS solution file ---*/
            SetCGNS_Solution(config[iZone], geometry[iZone], iZone);
            break;
            
          case PARAVIEW:
            
            /*--- Write a Paraview ASCII file ---*/
            SetParaview_ASCII(config[iZone], geometry[iZone], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone], false);
            break;
            
          default:
            break;
        }
        
      }
      
      if (Wrt_Srf) {
        
        if (rank == MASTER_NODE) cout <<"Writing surface solution file." << endl;
        
        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            SetTecplot_ASCII(config[iZone], geometry[iZone], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone], true);
            break;
            
          case TECPLOT_BINARY:
            
            /*--- Write a Tecplot binary solution file ---*/
            SetTecplot_SurfaceMesh(config[iZone], geometry[iZone], iZone);
            SetTecplot_SurfaceSolution(config[iZone], geometry[iZone], iZone);
            break;
            
          case PARAVIEW:
            
            /*--- Write a Paraview ASCII file ---*/
            SetParaview_ASCII(config[iZone], geometry[iZone], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone], true);
            break;
            
          default:
            break;
        }
      }
      
      if (FileFormat == TECPLOT_BINARY) {
        if (!wrote_base_file)
          DeallocateConnectivity(config[iZone], geometry[iZone], false);
        if (!wrote_surf_file)
          DeallocateConnectivity(config[iZone], geometry[iZone], wrote_surf_file);
      }
      
      if (Wrt_Vol || Wrt_Srf)
        DeallocateSolution(config[iZone], geometry[iZone]);
    }
    
    /*--- Final broadcast (informing other procs that the base output
     file was written) & barrier to sync up after master node writes
     output files. ---*/
    
#ifndef NO_MPI
    MPI::COMM_WORLD.Bcast(&wrote_base_file, 1, MPI::INT, MASTER_NODE);
    MPI::COMM_WORLD.Barrier();
#endif
    
  }
}

void COutput::SetEquivalentArea(CSolver *solver_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {
  
  ofstream EquivArea_file, FuncGrad_file;
  unsigned short iMarker = 0, iDim;
  short *AzimuthalAngle = NULL;
  double Gamma, auxXCoord, auxYCoord, auxZCoord, InverseDesign = 0.0, DeltaX, Coord_i, Coord_j, jp1Coord, *Coord = NULL, MeanFuntion,
  *Face_Normal = NULL, auxArea, auxPress, Mach, Beta, R_Plane, Pressure_Inf, Density_Inf,
  RefAreaCoeff, ModVelocity_Inf, Velocity_Inf[3], factor, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL,
  *Pressure = NULL, *FaceArea = NULL, *EquivArea = NULL, *TargetArea = NULL, *NearFieldWeight = NULL,
  *Weight = NULL, jFunction, jp1Function;
  unsigned long jVertex, iVertex, iPoint, nVertex_NearField = 0, auxPoint,
  *IdPoint = NULL, *IdDomain = NULL, auxDomain;
  unsigned short iPhiAngle;
  ofstream NearFieldEA_file; ifstream TargetEA_file;
  
  double XCoordBegin_OF = config->GetEA_IntLimit(0);
  double XCoordEnd_OF = config->GetEA_IntLimit(1);
  
  unsigned short nDim = geometry->GetnDim();
  double AoA = -(config->GetAoA()*PI_NUMBER/180.0);
  
  int rank = MESH_0;
  
  Mach  = config->GetMach_FreeStreamND();
  Gamma = config->GetGamma();
  Beta = sqrt(Mach*Mach-1.0);
  R_Plane = fabs(config->GetEA_IntLimit(2));
  Pressure_Inf = config->GetPressure_FreeStreamND();
  Density_Inf = config->GetDensity_FreeStreamND();
  RefAreaCoeff = config->GetRefAreaCoeff();
  Velocity_Inf[0] = config->GetVelocity_FreeStreamND()[0];
  Velocity_Inf[1] = config->GetVelocity_FreeStreamND()[1];
  Velocity_Inf[2] = config->GetVelocity_FreeStreamND()[2];
  ModVelocity_Inf = 0;
  for (iDim = 0; iDim < 3; iDim++)
    ModVelocity_Inf += Velocity_Inf[iDim] * Velocity_Inf[iDim];
  ModVelocity_Inf = sqrt(ModVelocity_Inf);
  
  factor = 4.0*sqrt(2.0*Beta*R_Plane) / (Gamma*Pressure_Inf*Mach*Mach);
  
#ifdef NO_MPI
  
  /*--- Compute the total number of points on the near-field ---*/
  nVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->node[iPoint]->GetCoord();
        /*--- Using Face_Normal(z), and Coord(z) we identify only a surface,
         note that there are 2 NEARFIELD_BOUNDARY surfaces ---*/
        if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) nVertex_NearField ++;
      }
  
  /*--- Create an array with all the coordinates, points, pressures, face area,
   equivalent area, and nearfield weight ---*/
  Xcoord = new double[nVertex_NearField];
  Ycoord = new double[nVertex_NearField];
  Zcoord = new double[nVertex_NearField];
  AzimuthalAngle = new short[nVertex_NearField];
  IdPoint = new unsigned long[nVertex_NearField];
  IdDomain = new unsigned long[nVertex_NearField];
  Pressure = new double[nVertex_NearField];
  FaceArea = new double[nVertex_NearField];
  EquivArea = new double[nVertex_NearField];
  TargetArea = new double[nVertex_NearField];
  NearFieldWeight = new double[nVertex_NearField];
  Weight = new double[nVertex_NearField];
  
  /*--- Copy the boundary information to an array ---*/
  nVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->node[iPoint]->GetCoord();
        
        if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) {
          
          IdPoint[nVertex_NearField] = iPoint;
          Xcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
          Ycoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(1);
          
          if (nDim ==2) {
            AzimuthalAngle[nVertex_NearField] = 0;
          }
          
          if (nDim == 3) {
            Zcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(2);
            
            /*--- Rotate the nearfield cylinder (AoA) only 3D ---*/
            double YcoordRot = Ycoord[nVertex_NearField];
            double ZcoordRot = Xcoord[nVertex_NearField]*sin(AoA) + Zcoord[nVertex_NearField]*cos(AoA);
            
            /* Compute the Azimuthal angle (resolution of degress in the Azimuthal angle)---*/
            double AngleDouble; short AngleInt;
            AngleDouble = atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER;
            AngleInt = (short) floor(AngleDouble + 0.5);
            if (AngleInt >= 0) AzimuthalAngle[nVertex_NearField] = AngleInt;
            else AzimuthalAngle[nVertex_NearField] = 180 + AngleInt;
          }
          
          if (AzimuthalAngle[nVertex_NearField] <= 60) {
            Pressure[nVertex_NearField] = solver_container->node[iPoint]->GetPressure(COMPRESSIBLE);
            FaceArea[nVertex_NearField] = fabs(Face_Normal[nDim-1]);
            nVertex_NearField ++;
          }
          
        }
      }
  
#else
  
  int nProcessor = MPI::COMM_WORLD.Get_size();
  rank = MPI::COMM_WORLD.Get_rank();
  
  unsigned long nLocalVertex_NearField = 0, MaxLocalVertex_NearField = 0;
  int iProcessor;
  
  unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];
  unsigned long *Buffer_Send_nVertex = new unsigned long [1];
  
  /*--- Compute the total number of points of the near-field ghost nodes ---*/
  nLocalVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->node[iPoint]->GetCoord();
        
        if (geometry->node[iPoint]->GetDomain())
          if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0))
            nLocalVertex_NearField ++;
      }
  
  Buffer_Send_nVertex[0] = nLocalVertex_NearField;
  
  /*--- Send Near-Field vertex information --*/
  MPI::COMM_WORLD.Allreduce(&nLocalVertex_NearField, &nVertex_NearField, 1, MPI::UNSIGNED_LONG, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI::UNSIGNED_LONG, MPI::MAX);
  MPI::COMM_WORLD.Allgather(Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG);
  
  double *Buffer_Send_Xcoord = new double[MaxLocalVertex_NearField];
  double *Buffer_Send_Ycoord = new double[MaxLocalVertex_NearField];
  double *Buffer_Send_Zcoord = new double[MaxLocalVertex_NearField];
  unsigned long *Buffer_Send_IdPoint = new unsigned long [MaxLocalVertex_NearField];
  double *Buffer_Send_Pressure = new double [MaxLocalVertex_NearField];
  double *Buffer_Send_FaceArea = new double[MaxLocalVertex_NearField];
  
  double *Buffer_Receive_Xcoord = new double[nProcessor*MaxLocalVertex_NearField];
  double *Buffer_Receive_Ycoord = new double[nProcessor*MaxLocalVertex_NearField];
  double *Buffer_Receive_Zcoord = new double[nProcessor*MaxLocalVertex_NearField];
  unsigned long *Buffer_Receive_IdPoint = new unsigned long[nProcessor*MaxLocalVertex_NearField];
  double *Buffer_Receive_Pressure = new double[nProcessor*MaxLocalVertex_NearField];
  double *Buffer_Receive_FaceArea = new double[nProcessor*MaxLocalVertex_NearField];
  
  unsigned long nBuffer_Xcoord = MaxLocalVertex_NearField;
  unsigned long nBuffer_Ycoord = MaxLocalVertex_NearField;
  unsigned long nBuffer_Zcoord = MaxLocalVertex_NearField;
  unsigned long nBuffer_IdPoint = MaxLocalVertex_NearField;
  unsigned long nBuffer_Pressure = MaxLocalVertex_NearField;
  unsigned long nBuffer_FaceArea = MaxLocalVertex_NearField;
  
  for (iVertex = 0; iVertex < MaxLocalVertex_NearField; iVertex++) {
    Buffer_Send_IdPoint[iVertex] = 0; Buffer_Send_Pressure[iVertex] = 0.0;
    Buffer_Send_FaceArea[iVertex] = 0.0; Buffer_Send_Xcoord[iVertex] = 0.0;
    Buffer_Send_Ycoord[iVertex] = 0.0; Buffer_Send_Zcoord[iVertex] = 0.0;
  }
  
  /*--- Copy coordinates, index points, and pressures to the auxiliar vector --*/
  nLocalVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->node[iPoint]->GetCoord();
        
        if (geometry->node[iPoint]->GetDomain())
          if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) {
            Buffer_Send_IdPoint[nLocalVertex_NearField] = iPoint;
            Buffer_Send_Xcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
            Buffer_Send_Ycoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(1);
            Buffer_Send_Zcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(2);
            Buffer_Send_Pressure[nLocalVertex_NearField] = solver_container->node[iPoint]->GetPressure(COMPRESSIBLE);
            Buffer_Send_FaceArea[nLocalVertex_NearField] = fabs(Face_Normal[nDim-1]);
            nLocalVertex_NearField++;
          }
      }
  
  /*--- Send all the information --*/
  MPI::COMM_WORLD.Gather(Buffer_Send_Xcoord, nBuffer_Xcoord, MPI::DOUBLE, Buffer_Receive_Xcoord, nBuffer_Xcoord, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_Ycoord, nBuffer_Ycoord, MPI::DOUBLE, Buffer_Receive_Ycoord, nBuffer_Ycoord, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_Zcoord, nBuffer_Zcoord, MPI::DOUBLE, Buffer_Receive_Zcoord, nBuffer_Zcoord, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_IdPoint, nBuffer_IdPoint, MPI::UNSIGNED_LONG, Buffer_Receive_IdPoint, nBuffer_IdPoint, MPI::UNSIGNED_LONG, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_Pressure, nBuffer_Pressure, MPI::DOUBLE, Buffer_Receive_Pressure, nBuffer_Pressure, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_FaceArea, nBuffer_FaceArea, MPI::DOUBLE, Buffer_Receive_FaceArea, nBuffer_FaceArea, MPI::DOUBLE, MASTER_NODE);
  
  if (rank == MASTER_NODE) {
    
    Xcoord = new double[nVertex_NearField];
    Ycoord = new double[nVertex_NearField];
    Zcoord = new double[nVertex_NearField];
    AzimuthalAngle = new short[nVertex_NearField];
    IdPoint = new unsigned long[nVertex_NearField];
    IdDomain = new unsigned long[nVertex_NearField];
    Pressure = new double[nVertex_NearField];
    FaceArea = new double[nVertex_NearField];
    EquivArea = new double[nVertex_NearField];
    TargetArea = new double[nVertex_NearField];
    NearFieldWeight = new double[nVertex_NearField];
    Weight = new double[nVertex_NearField];
    
    nVertex_NearField = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
      for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
        Xcoord[nVertex_NearField] = Buffer_Receive_Xcoord[iProcessor*MaxLocalVertex_NearField+iVertex];
        Ycoord[nVertex_NearField] = Buffer_Receive_Ycoord[iProcessor*MaxLocalVertex_NearField+iVertex];
        
        if (nDim == 2) {
          AzimuthalAngle[nVertex_NearField] = 0;
        }
        
        if (nDim == 3) {
          Zcoord[nVertex_NearField] = Buffer_Receive_Zcoord[iProcessor*MaxLocalVertex_NearField+iVertex];
          
          /*--- Rotate the nearfield cylinder  ---*/
          double YcoordRot = Ycoord[nVertex_NearField];
          double ZcoordRot = Xcoord[nVertex_NearField]*sin(AoA) + Zcoord[nVertex_NearField]*cos(AoA);
          
          /*--- Compute the Azimuthal angle ---*/
          double AngleDouble; short AngleInt;
          AngleDouble = atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER;
          AngleInt = (short) floor(AngleDouble + 0.5);
          if (AngleInt >= 0) AzimuthalAngle[nVertex_NearField] = AngleInt;
          else AzimuthalAngle[nVertex_NearField] = 180 + AngleInt;
        }
        
        if (AzimuthalAngle[nVertex_NearField] <= 60) {
          IdPoint[nVertex_NearField] = Buffer_Receive_IdPoint[iProcessor*MaxLocalVertex_NearField+iVertex];
          Pressure[nVertex_NearField] = Buffer_Receive_Pressure[iProcessor*MaxLocalVertex_NearField+iVertex];
          FaceArea[nVertex_NearField] = Buffer_Receive_FaceArea[iProcessor*MaxLocalVertex_NearField+iVertex];
          IdDomain[nVertex_NearField] = iProcessor;
          nVertex_NearField++;
        }
        
      }
  }
  
  delete [] Buffer_Receive_nVertex;
  delete [] Buffer_Send_nVertex;
  
  delete [] Buffer_Send_Xcoord;
  delete [] Buffer_Send_Ycoord;
  delete [] Buffer_Send_Zcoord;
  delete [] Buffer_Send_IdPoint;
  delete [] Buffer_Send_Pressure;
  delete [] Buffer_Send_FaceArea;
  
  delete [] Buffer_Receive_Xcoord;
  delete [] Buffer_Receive_IdPoint;
  delete [] Buffer_Receive_Pressure;
  delete [] Buffer_Receive_FaceArea;
  
#endif
  
  if (rank == MASTER_NODE) {
    
    vector<short> PhiAngleList;
    vector<short>::iterator IterPhiAngleList;
    
    for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
      PhiAngleList.push_back(AzimuthalAngle[iVertex]);
    
    sort( PhiAngleList.begin(), PhiAngleList.end());
    IterPhiAngleList = unique( PhiAngleList.begin(), PhiAngleList.end());
    PhiAngleList.resize( IterPhiAngleList - PhiAngleList.begin() );
    
    /*--- Create vectors and distribute the values among the different PhiAngle queues ---*/
    vector<vector<double> > Xcoord_PhiAngle; Xcoord_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > Ycoord_PhiAngle; Ycoord_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > Zcoord_PhiAngle; Zcoord_PhiAngle.resize(PhiAngleList.size());
    vector<vector<unsigned long> > IdPoint_PhiAngle; IdPoint_PhiAngle.resize(PhiAngleList.size());
    vector<vector<unsigned long> > IdDomain_PhiAngle; IdDomain_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > Pressure_PhiAngle; Pressure_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > FaceArea_PhiAngle; FaceArea_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > EquivArea_PhiAngle; EquivArea_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > TargetArea_PhiAngle; TargetArea_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > NearFieldWeight_PhiAngle; NearFieldWeight_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > Weight_PhiAngle; Weight_PhiAngle.resize(PhiAngleList.size());
    
    /*--- Distribute the values among the different PhiAngles ---*/
    for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        if (AzimuthalAngle[iVertex] == PhiAngleList[iPhiAngle]) {
          Xcoord_PhiAngle[iPhiAngle].push_back(Xcoord[iVertex]);
          Ycoord_PhiAngle[iPhiAngle].push_back(Ycoord[iVertex]);
          Zcoord_PhiAngle[iPhiAngle].push_back(Zcoord[iVertex]);
          IdPoint_PhiAngle[iPhiAngle].push_back(IdPoint[iVertex]);
          IdDomain_PhiAngle[iPhiAngle].push_back(IdDomain[iVertex]);
          Pressure_PhiAngle[iPhiAngle].push_back(Pressure[iVertex]);
          FaceArea_PhiAngle[iPhiAngle].push_back(FaceArea[iVertex]);
          EquivArea_PhiAngle[iPhiAngle].push_back(EquivArea[iVertex]);
          TargetArea_PhiAngle[iPhiAngle].push_back(TargetArea[iVertex]);
          NearFieldWeight_PhiAngle[iPhiAngle].push_back(NearFieldWeight[iVertex]);
          Weight_PhiAngle[iPhiAngle].push_back(Weight[iVertex]);
        }
    
    /*--- Order the arrays (x Coordinate, Pressure, Point, and Domain) ---*/
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      for (iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++)
        for (jVertex = 0; jVertex < Xcoord_PhiAngle[iPhiAngle].size() - 1 - iVertex; jVertex++)
          if (Xcoord_PhiAngle[iPhiAngle][jVertex] > Xcoord_PhiAngle[iPhiAngle][jVertex+1]) {
            auxXCoord = Xcoord_PhiAngle[iPhiAngle][jVertex]; Xcoord_PhiAngle[iPhiAngle][jVertex] = Xcoord_PhiAngle[iPhiAngle][jVertex+1]; Xcoord_PhiAngle[iPhiAngle][jVertex+1] = auxXCoord;
            auxYCoord = Ycoord_PhiAngle[iPhiAngle][jVertex]; Ycoord_PhiAngle[iPhiAngle][jVertex] = Ycoord_PhiAngle[iPhiAngle][jVertex+1]; Ycoord_PhiAngle[iPhiAngle][jVertex+1] = auxYCoord;
            auxZCoord = Zcoord_PhiAngle[iPhiAngle][jVertex]; Zcoord_PhiAngle[iPhiAngle][jVertex] = Zcoord_PhiAngle[iPhiAngle][jVertex+1]; Zcoord_PhiAngle[iPhiAngle][jVertex+1] = auxZCoord;
            auxPress = Pressure_PhiAngle[iPhiAngle][jVertex]; Pressure_PhiAngle[iPhiAngle][jVertex] = Pressure_PhiAngle[iPhiAngle][jVertex+1]; Pressure_PhiAngle[iPhiAngle][jVertex+1] = auxPress;
            auxArea = FaceArea_PhiAngle[iPhiAngle][jVertex]; FaceArea_PhiAngle[iPhiAngle][jVertex] = FaceArea_PhiAngle[iPhiAngle][jVertex+1]; FaceArea_PhiAngle[iPhiAngle][jVertex+1] = auxArea;
            auxPoint = IdPoint_PhiAngle[iPhiAngle][jVertex]; IdPoint_PhiAngle[iPhiAngle][jVertex] = IdPoint_PhiAngle[iPhiAngle][jVertex+1]; IdPoint_PhiAngle[iPhiAngle][jVertex+1] = auxPoint;
            auxDomain = IdDomain_PhiAngle[iPhiAngle][jVertex]; IdDomain_PhiAngle[iPhiAngle][jVertex] = IdDomain_PhiAngle[iPhiAngle][jVertex+1]; IdDomain_PhiAngle[iPhiAngle][jVertex+1] = auxDomain;
          }
    
    
    /*--- Check that all the azimuth lists have the same size ---*/
    unsigned short nVertex = Xcoord_PhiAngle[0].size();
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      unsigned short nVertex_aux = Xcoord_PhiAngle[iPhiAngle].size();
      if (nVertex_aux != nVertex) cout <<"Be careful!!! one azimuth list is shorter than the other"<< endl;
      nVertex = min(nVertex, nVertex_aux);
    }
    
    /*--- Compute equivalent area distribution at each azimuth angle ---*/
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      EquivArea_PhiAngle[iPhiAngle][0] = 0.0;
      for (iVertex = 1; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
        EquivArea_PhiAngle[iPhiAngle][iVertex] = 0.0;
        
        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][iVertex]*sin(AoA);
        
        for (jVertex = 0; jVertex < iVertex-1; jVertex++) {
          
          Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][jVertex]*sin(AoA);
          jp1Coord = Xcoord_PhiAngle[iPhiAngle][jVertex+1]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][jVertex+1]*sin(AoA);
          
          jFunction = factor*(Pressure_PhiAngle[iPhiAngle][jVertex] - Pressure_Inf)*sqrt(Coord_i-Coord_j);
          jp1Function = factor*(Pressure_PhiAngle[iPhiAngle][jVertex+1] - Pressure_Inf)*sqrt(Coord_i-jp1Coord);
          
          DeltaX = (jp1Coord-Coord_j);
          MeanFuntion = 0.5*(jp1Function + jFunction);
          EquivArea_PhiAngle[iPhiAngle][iVertex] += DeltaX * MeanFuntion;
        }
      }
    }
    
    /*--- Create a file with the equivalent area distribution at each azimuthal angle ---*/
    NearFieldEA_file.precision(15);
    NearFieldEA_file.open("NearFieldEA.plt", ios::out);
    NearFieldEA_file << "TITLE = \"Nearfield Equivalent Area at each azimuthal angle \"" << endl;
    NearFieldEA_file << "VARIABLES = \"Coord (local to the near-field cylinder)\"";
    
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      NearFieldEA_file << ", \"Equivalent Area, Phi= " << PhiAngleList[iPhiAngle] << " deg.\"";
    }
    
    NearFieldEA_file << endl;
    for (iVertex = 0; iVertex < EquivArea_PhiAngle[0].size(); iVertex++) {
      double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
      double XcoordRot_init = Xcoord_PhiAngle[0][0]*cos(AoA) - Zcoord_PhiAngle[0][0]*sin(AoA);
      NearFieldEA_file << scientific << XcoordRot - XcoordRot_init;
      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
        NearFieldEA_file << scientific << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex];
      }
      NearFieldEA_file << endl;
    }
    NearFieldEA_file.close();
    
    /*--- Read target equivalent area from the configuration file,
     this first implementation requires a complete table (same as the original
     EA table). so... no interpolation. ---*/
    
    vector<vector<double> > TargetArea_PhiAngle_Trans;
    TargetEA_file.open("TargetEA.dat", ios::in);
    
    if (TargetEA_file.fail()) {
      if (iExtIter == 0) { cout << "There is no Target Equivalent Area file (TargetEA.dat)!!"<< endl;
        cout << "Using default parameters (Target Equiv Area = 0.0)" << endl;
      }
      /*--- Set the table to 0 ---*/
      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        for (iVertex = 0; iVertex < TargetArea_PhiAngle[iPhiAngle].size(); iVertex++)
          TargetArea_PhiAngle[iPhiAngle][iVertex] = 0.0;
    }
    else {
      
      /*--- skip header lines ---*/
      string line;
      getline(TargetEA_file, line);
      getline(TargetEA_file, line);
      
      while (TargetEA_file) {
        
        string line;
        getline(TargetEA_file, line);
        istringstream is(line);
        vector<double> row;
        unsigned short iter = 0;
        
        while (is.good()) {
          string token;
          getline(is,token,',');
          
          istringstream js(token);
          
          double data;
          js >> data;
          
          /*--- The first element in the table is the coordinate ---*/
          if (iter != 0) row.push_back(data);
          iter++;
        }
        TargetArea_PhiAngle_Trans.push_back(row);
      }
      
      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++)
          TargetArea_PhiAngle[iPhiAngle][iVertex] = TargetArea_PhiAngle_Trans[iVertex][iPhiAngle];
    }
    
    /*--- Divide by the number of Phi angles in the nearfield ---*/
    double PhiFactor = 1.0/double(PhiAngleList.size());
    
    /*--- Evaluate the objective function ---*/
    InverseDesign = 0;
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
        Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;
        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];
        
        double Difference = EquivArea_PhiAngle[iPhiAngle][iVertex]-TargetArea_PhiAngle[iPhiAngle][iVertex];
        if ((Coord_i < XCoordBegin_OF) || (Coord_i > XCoordEnd_OF)) Difference = 0.0;
        
        InverseDesign += PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex]*Difference*Difference;
        
      }
    
    /*--- Evaluate the weight of the nearfield pressure (adjoint input) ---*/
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];
        NearFieldWeight_PhiAngle[iPhiAngle][iVertex] = 0.0;
        for (jVertex = iVertex; jVertex < EquivArea_PhiAngle[iPhiAngle].size(); jVertex++) {
          Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex];
          Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;
          
          double Difference = EquivArea_PhiAngle[iPhiAngle][jVertex]-TargetArea_PhiAngle[iPhiAngle][jVertex];
          if ((Coord_j < XCoordBegin_OF) || (Coord_j > XCoordEnd_OF)) Difference = 0.0;
          
          NearFieldWeight_PhiAngle[iPhiAngle][iVertex] += PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex]*2.0*Difference*factor*sqrt(Coord_j-Coord_i);
        }
      }		
    
    /*--- Write the Nearfield pressure at each Azimuthal PhiAngle ---*/
    EquivArea_file.precision(15);
    EquivArea_file.open("nearfield_flow.plt", ios::out);
    EquivArea_file << "TITLE = \"SU2 Equivalent area computation at each azimuthal angle \"" << endl;
    EquivArea_file << "VARIABLES = \"Coord (local to the near-field cylinder)\",\"Equivalent area\",\"Target equivalent area\",\"NearField weight\",\"Pressure coefficient\"" << endl;
    
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      EquivArea_file << fixed << "ZONE T= \"Azimuthal angle " << PhiAngleList[iPhiAngle] << " deg.\"" << endl;
      for (iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++) {
        
        double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
        double XcoordRot_init = Xcoord_PhiAngle[0][0]*cos(AoA) - Zcoord_PhiAngle[0][0]*sin(AoA);
        
        EquivArea_file << scientific << XcoordRot - XcoordRot_init << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex]
        << ", " << TargetArea_PhiAngle[iPhiAngle][iVertex] << ", " << NearFieldWeight_PhiAngle[iPhiAngle][iVertex] << ", " <<
        (Pressure_PhiAngle[iPhiAngle][iVertex]-Pressure_Inf)/Pressure_Inf << endl;
      }
    }
    
    EquivArea_file.close();
    
    /*--- Write Weight file for adjoint computation ---*/
    FuncGrad_file.precision(15);
    FuncGrad_file.open("WeightNF.dat", ios::out);
    
    FuncGrad_file << scientific << "-1.0";
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      FuncGrad_file << scientific << "\t" << PhiAngleList[iPhiAngle];
    FuncGrad_file << endl;
    
    for (iVertex = 0; iVertex < NearFieldWeight_PhiAngle[0].size(); iVertex++) {
      double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
      FuncGrad_file << scientific << XcoordRot;
      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        FuncGrad_file << scientific << "\t" << NearFieldWeight_PhiAngle[iPhiAngle][iVertex];
      FuncGrad_file << endl;
    }
    FuncGrad_file.close();		
    
    /*--- Delete structures ---*/
    delete [] Xcoord; delete [] Ycoord; delete [] Zcoord; 
    delete [] AzimuthalAngle; delete [] IdPoint; delete [] IdDomain;
    delete [] Pressure; delete [] FaceArea;
    delete [] EquivArea; delete [] TargetArea;
    delete [] NearFieldWeight; delete [] Weight;
    
  }
  
#ifdef NO_MPI
  
  /*--- Store the value of the NearField coefficient ---*/
  solver_container->SetTotal_CEquivArea(InverseDesign);
  
#else
  
  /*--- Send the value of the NearField coefficient to all the processors ---*/
  MPI::COMM_WORLD.Bcast (&InverseDesign, 1, MPI::DOUBLE, MASTER_NODE);
  
  /*--- Store the value of the NearField coefficient ---*/
  solver_container->SetTotal_CEquivArea(InverseDesign);
  
#endif
  
}
