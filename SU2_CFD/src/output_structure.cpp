/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information
 * \author F. Palacios, T. Economon
 * \version 4.0.2 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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
  nGlobal_Pris      = 0;
  nGlobal_Pyra      = 0;
  nGlobal_Line      = 0;
  nGlobal_BoundTria = 0;
  nGlobal_BoundQuad = 0;
  
  /*--- Initialize CGNS write flag ---*/
  
  wrote_base_file = false;
  
  /*--- Initialize CGNS write flag ---*/
  
  wrote_CGNS_base = false;
  
  /*--- Initialize Tecplot surface flag ---*/
  
  wrote_surf_file = false;
  
  /*--- Initialize Paraview write flag ---*/
  
  wrote_Paraview_base = false;
  
  /*--- Initialize residual ---*/

  RhoRes_New = EPS;
  RhoRes_Old = EPS;
  
}

COutput::~COutput(void) { }

void COutput::SetSurfaceCSV_Flow(CConfig *config, CGeometry *geometry,
                                 CSolver *FlowSolver, unsigned long iExtIter,
                                 unsigned short val_iZone) {
  
  unsigned short iMarker;
  unsigned long iPoint, iVertex, Global_Index;
  su2double PressCoeff = 0.0, SkinFrictionCoeff;
  su2double xCoord = 0.0, yCoord = 0.0, zCoord = 0.0, Mach, Pressure;
  char cstr[200];
  
  unsigned short solver = config->GetKind_Solver();
  unsigned short nDim = geometry->GetnDim();
  
#ifndef HAVE_MPI
  
  su2double HeatFlux;
  char buffer [50];
  ofstream SurfFlow_file;
  
  /*--- Write file name with extension if unsteady ---*/
  strcpy (cstr, config->GetSurfFlowCoeff_FileName().c_str());
  
  if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    if (SU2_TYPE::Int(val_iZone) < 10) SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(val_iZone));
    if ((SU2_TYPE::Int(val_iZone) >= 10)   && (SU2_TYPE::Int(val_iZone) < 100))   SPRINTF (buffer, "_000%d.csv", SU2_TYPE::Int(val_iZone));
    if ((SU2_TYPE::Int(val_iZone) >= 100)  && (SU2_TYPE::Int(val_iZone) < 1000))  SPRINTF (buffer, "_00%d.csv", SU2_TYPE::Int(val_iZone));
    if ((SU2_TYPE::Int(val_iZone) >= 1000) && (SU2_TYPE::Int(val_iZone) < 10000)) SPRINTF (buffer, "_0%d.csv", SU2_TYPE::Int(val_iZone));
    if (SU2_TYPE::Int(val_iZone) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(val_iZone));
    
  } else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.csv",  SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.csv",   SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.csv",    SU2_TYPE::Int(iExtIter));
    if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(iExtIter));
  }
  else
    SPRINTF (buffer, ".csv");
  
  strcat (cstr, buffer);
  SurfFlow_file.precision(15);
  SurfFlow_file.open(cstr, ios::out);
  
  SurfFlow_file << "\"Global_Index\", \"x_coord\", \"y_coord\", ";
  if (nDim == 3) SurfFlow_file << "\"z_coord\", ";
  SurfFlow_file << "\"Pressure\", \"Pressure_Coefficient\", ";
  
  switch (solver) {
    case EULER : SurfFlow_file <<  "\"Mach_Number\"" << endl; break;
    case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"Skin_Friction_Coefficient\", \"Heat_Flux\"" << endl; break;
  }
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        xCoord = geometry->node[iPoint]->GetCoord(0);
        yCoord = geometry->node[iPoint]->GetCoord(1);
        if (nDim == 3) zCoord = geometry->node[iPoint]->GetCoord(2);
        
        /*--- The output should be in inches ---*/
        
        if (config->GetSystemMeasurements() == US) {
          xCoord *= 12.0; yCoord *= 12.0;
          if (nDim == 3) zCoord *= 12.0;
        }
        
        Pressure = FlowSolver->node[iPoint]->GetPressure();
        PressCoeff = FlowSolver->GetCPressure(iMarker, iVertex);
        SurfFlow_file << scientific << Global_Index << ", " << xCoord << ", " << yCoord << ", ";
        if (nDim == 3) SurfFlow_file << scientific << zCoord << ", ";
        SurfFlow_file << scientific << Pressure << ", " << PressCoeff << ", ";
        switch (solver) {
          case EULER :
            Mach = sqrt(FlowSolver->node[iPoint]->GetVelocity2()) / FlowSolver->node[iPoint]->GetSoundSpeed();
            SurfFlow_file << scientific << Mach << endl;
            break;
          case NAVIER_STOKES: case RANS:
            SkinFrictionCoeff = FlowSolver->GetCSkinFriction(iMarker, iVertex);
            HeatFlux = FlowSolver->GetHeatFlux(iMarker, iVertex);
            SurfFlow_file << scientific << SkinFrictionCoeff << ", " << HeatFlux << endl;
            break;
        }
      }
    }
  }
  
  SurfFlow_file.close();
  
#else
  
  int rank, iProcessor, nProcessor;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  unsigned long Buffer_Send_nVertex[1], *Buffer_Recv_nVertex = NULL;
  unsigned long nVertex_Surface = 0, nLocalVertex_Surface = 0;
  unsigned long MaxLocalVertex_Surface = 0;
  
  /*--- Find the max number of surface vertices among all
   partitions and set up buffers. The master node will handle the
   writing of the CSV file after gathering all of the data. ---*/
  
  nLocalVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface++;
      }
  
  /*--- Communicate the number of local vertices on each partition
   to the master node ---*/
  
  Buffer_Send_nVertex[0] = nLocalVertex_Surface;
  if (rank == MASTER_NODE) Buffer_Recv_nVertex = new unsigned long [nProcessor];
  
  SU2_MPI::Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- Send and Recv buffers ---*/
  
  su2double *Buffer_Send_Coord_x = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_Coord_x = NULL;
  
  su2double *Buffer_Send_Coord_y = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_Coord_y = NULL;
  
  su2double *Buffer_Send_Coord_z = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_Coord_z = NULL;
  
  su2double *Buffer_Send_Press = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_Press = NULL;
  
  su2double *Buffer_Send_CPress = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_CPress = NULL;
  
  su2double *Buffer_Send_Mach = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_Mach = NULL;
  
  su2double *Buffer_Send_SkinFriction = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_SkinFriction = NULL;
  
  su2double *Buffer_Send_HeatTransfer = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_HeatTransfer = NULL;
  
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long [MaxLocalVertex_Surface];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Prepare the receive buffers on the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    Buffer_Recv_Coord_x = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_Coord_y = new su2double [nProcessor*MaxLocalVertex_Surface];
    if (nDim == 3) Buffer_Recv_Coord_z = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_Press   = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_CPress  = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_Mach    = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_SkinFriction = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_HeatTransfer = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_GlobalIndex  = new unsigned long [nProcessor*MaxLocalVertex_Surface];
  }
  
  /*--- Loop over all vertices in this partition and load the
   data of the specified type into the buffer to be sent to
   the master node. ---*/
  
  nVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          Buffer_Send_Press[nVertex_Surface] = FlowSolver->node[iPoint]->GetPressure();
          Buffer_Send_CPress[nVertex_Surface] = FlowSolver->GetCPressure(iMarker, iVertex);
          Buffer_Send_Coord_x[nVertex_Surface] = geometry->node[iPoint]->GetCoord(0);
          Buffer_Send_Coord_y[nVertex_Surface] = geometry->node[iPoint]->GetCoord(1);
          if (nDim == 3) { Buffer_Send_Coord_z[nVertex_Surface] = geometry->node[iPoint]->GetCoord(2); }
          
          /*--- If US system, the output should be in inches ---*/
          
          if (config->GetSystemMeasurements() == US) {
            Buffer_Send_Coord_x[nVertex_Surface] *= 12.0;
            Buffer_Send_Coord_y[nVertex_Surface] *= 12.0;
            if (nDim == 3) Buffer_Send_Coord_z[nVertex_Surface] *= 12.0;
          }
          
          Buffer_Send_GlobalIndex[nVertex_Surface] = geometry->node[iPoint]->GetGlobalIndex();
          
          if (solver == EULER)
            Buffer_Send_Mach[nVertex_Surface] = sqrt(FlowSolver->node[iPoint]->GetVelocity2()) / FlowSolver->node[iPoint]->GetSoundSpeed();
          if ((solver == NAVIER_STOKES) || (solver == RANS))
            Buffer_Send_SkinFriction[nVertex_Surface] = FlowSolver->GetCSkinFriction(iMarker, iVertex);
          nVertex_Surface++;
        }
      }
  
  /*--- Send the information to the master node ---*/
  
  SU2_MPI::Gather(Buffer_Send_Coord_x, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_x, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Coord_y, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_y, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) SU2_MPI::Gather(Buffer_Send_Coord_z, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_z, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Press, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Press, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_CPress, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_CPress, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (solver == EULER) SU2_MPI::Gather(Buffer_Send_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if ((solver == NAVIER_STOKES) || (solver == RANS)) SU2_MPI::Gather(Buffer_Send_SkinFriction, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_SkinFriction, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_GlobalIndex, MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- The master node unpacks the data and writes the surface CSV file ---*/
  
  if (rank == MASTER_NODE) {
    
    /*--- Write file name with extension if unsteady ---*/
    char buffer[50];
    string filename = config->GetSurfFlowCoeff_FileName();
    ofstream SurfFlow_file;
    
    /*--- Write file name with extension if unsteady ---*/
    strcpy (cstr, filename.c_str());
    if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
      if (SU2_TYPE::Int(val_iZone) < 10) SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(val_iZone));
      if ((SU2_TYPE::Int(val_iZone) >= 10) && (SU2_TYPE::Int(val_iZone) < 100)) SPRINTF (buffer, "_000%d.csv", SU2_TYPE::Int(val_iZone));
      if ((SU2_TYPE::Int(val_iZone) >= 100) && (SU2_TYPE::Int(val_iZone) < 1000)) SPRINTF (buffer, "_00%d.csv", SU2_TYPE::Int(val_iZone));
      if ((SU2_TYPE::Int(val_iZone) >= 1000) && (SU2_TYPE::Int(val_iZone) < 10000)) SPRINTF (buffer, "_0%d.csv", SU2_TYPE::Int(val_iZone));
      if (SU2_TYPE::Int(val_iZone) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(val_iZone));
      
    } else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
      if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.csv",  SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.csv",   SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.csv",    SU2_TYPE::Int(iExtIter));
      if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(iExtIter));
    }
    else
      SPRINTF (buffer, ".csv");
    
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
    
    /*--- Loop through all of the collected data and write each node's values ---*/
    
    unsigned long Total_Index;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {
        
        /*--- Current index position and global index ---*/
        Total_Index  = iProcessor*MaxLocalVertex_Surface+iVertex;
        Global_Index = Buffer_Recv_GlobalIndex[Total_Index];
        
        /*--- Retrieve the merged data for this node ---*/
        xCoord = Buffer_Recv_Coord_x[Total_Index];
        yCoord = Buffer_Recv_Coord_y[Total_Index];
        if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];
        Pressure   = Buffer_Recv_Press[Total_Index];
        PressCoeff = Buffer_Recv_CPress[Total_Index];
        
        /*--- Write the first part of the data ---*/
        SurfFlow_file << scientific << Global_Index << ", " << xCoord << ", " << yCoord << ", ";
        if (nDim == 3) SurfFlow_file << scientific << zCoord << ", ";
        SurfFlow_file << scientific << Pressure << ", " << PressCoeff << ", ";
        
        /*--- Write the solver-dependent part of the data ---*/
        switch (solver) {
          case EULER :
            Mach = Buffer_Recv_Mach[Total_Index];
            SurfFlow_file << scientific << Mach << endl;
            break;
          case NAVIER_STOKES: case RANS:
            SkinFrictionCoeff = Buffer_Recv_SkinFriction[Total_Index];
            SurfFlow_file << scientific << SkinFrictionCoeff << endl;
            break;
        }
      }
    }
    
    /*--- Close the CSV file ---*/
    SurfFlow_file.close();
    
    /*--- Release the recv buffers on the master node ---*/
    
    delete [] Buffer_Recv_Coord_x;
    delete [] Buffer_Recv_Coord_y;
    if (nDim == 3) delete [] Buffer_Recv_Coord_z;
    delete [] Buffer_Recv_Press;
    delete [] Buffer_Recv_CPress;
    delete [] Buffer_Recv_Mach;
    delete [] Buffer_Recv_SkinFriction;
    delete [] Buffer_Recv_HeatTransfer;
    delete [] Buffer_Recv_GlobalIndex;
    
    delete [] Buffer_Recv_nVertex;
    
  }
  
  /*--- Release the memory for the remaining buffers and exit ---*/
  
  delete [] Buffer_Send_Coord_x;
  delete [] Buffer_Send_Coord_y;
  delete [] Buffer_Send_Coord_z;
  delete [] Buffer_Send_Press;
  delete [] Buffer_Send_CPress;
  delete [] Buffer_Send_Mach;
  delete [] Buffer_Send_SkinFriction;
  delete [] Buffer_Send_HeatTransfer;
  delete [] Buffer_Send_GlobalIndex;
  
#endif
  
}

void COutput::SetSurfaceCSV_Adjoint(CConfig *config, CGeometry *geometry, CSolver *AdjSolver, CSolver *FlowSolution, unsigned long iExtIter, unsigned short val_iZone) {
  
#ifndef HAVE_MPI
  
  unsigned long iPoint, iVertex, Global_Index;
  su2double *Solution, xCoord, yCoord, zCoord;
  unsigned short iMarker;
  char cstr[200], buffer[50];
  ofstream SurfAdj_file;
  
  /*--- Write file name with extension if unsteady ---*/
  strcpy (cstr, config->GetSurfAdjCoeff_FileName().c_str());
  
  if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    if (SU2_TYPE::Int(val_iZone) < 10) SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(val_iZone));
    if ((SU2_TYPE::Int(val_iZone) >= 10) && (SU2_TYPE::Int(val_iZone) < 100)) SPRINTF (buffer, "_000%d.csv", SU2_TYPE::Int(val_iZone));
    if ((SU2_TYPE::Int(val_iZone) >= 100) && (SU2_TYPE::Int(val_iZone) < 1000)) SPRINTF (buffer, "_00%d.csv", SU2_TYPE::Int(val_iZone));
    if ((SU2_TYPE::Int(val_iZone) >= 1000) && (SU2_TYPE::Int(val_iZone) < 10000)) SPRINTF (buffer, "_0%d.csv", SU2_TYPE::Int(val_iZone));
    if (SU2_TYPE::Int(val_iZone) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(val_iZone));
    
  } else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.csv",  SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.csv",   SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.csv",    SU2_TYPE::Int(iExtIter));
    if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(iExtIter));
  }
  else
    SPRINTF (buffer, ".csv");
  
  strcat(cstr, buffer);
  SurfAdj_file.precision(15);
  SurfAdj_file.open(cstr, ios::out);
  
  if (geometry->GetnDim() == 2) {
    SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"";
    if (config->GetDiscrete_Adjoint()){
      SurfAdj_file << ",\"x_Sens\",\"y_Sens\"";
    }
    SurfAdj_file << endl;

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Global_Index = geometry->node[iPoint]->GetGlobalIndex();
          Solution = AdjSolver->node[iPoint]->GetSolution();
          xCoord = geometry->node[iPoint]->GetCoord(0);
          yCoord = geometry->node[iPoint]->GetCoord(1);
          
          /*--- If US system, the output should be in inches ---*/
          
          if (config->GetSystemMeasurements() == US) {
            xCoord *= 12.0;
            yCoord *= 12.0;
          }
          
          SurfAdj_file << scientific << Global_Index << ", " << AdjSolver->GetCSensitivity(iMarker, iVertex) << ", " << Solution[0] << ", "
          << Solution[1] << ", " << Solution[2] << ", " << Solution[3] <<", " << xCoord <<", "<< yCoord;
          if (config->GetDiscrete_Adjoint()){
            SurfAdj_file << ", " << AdjSolver->node[iPoint]->GetSensitivity(0) << ", " << AdjSolver->node[iPoint]->GetSensitivity(1);
          }
          SurfAdj_file << endl;
        }
    }
  }
  
  if (geometry->GetnDim() == 3) {
    SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"";
    if (config->GetDiscrete_Adjoint()){
      SurfAdj_file << ",\"x_Sens\",\"y_Sens\",\"z_Sens\"";
    }
    SurfAdj_file << endl;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Global_Index = geometry->node[iPoint]->GetGlobalIndex();
          Solution = AdjSolver->node[iPoint]->GetSolution();
          
          xCoord = geometry->node[iPoint]->GetCoord(0);
          yCoord = geometry->node[iPoint]->GetCoord(1);
          zCoord = geometry->node[iPoint]->GetCoord(2);
          
          /*--- If US system, the output should be in inches ---*/
          
          if (config->GetSystemMeasurements() == US) {
            xCoord *= 12.0;
            yCoord *= 12.0;
            zCoord *= 12.0;
          }
          
          SurfAdj_file << scientific << Global_Index << ", " << AdjSolver->GetCSensitivity(iMarker, iVertex) << ", " << Solution[0] << ", "
          << Solution[1] << ", " << Solution[2] << ", " << Solution[3] << ", " << Solution[4] << ", "<< xCoord <<", "<< yCoord <<", "<< zCoord;
          if (config->GetDiscrete_Adjoint()){
            SurfAdj_file << ", " << AdjSolver->node[iPoint]->GetSensitivity(0) << ", " << AdjSolver->node[iPoint]->GetSensitivity(1)
                         << ", " << AdjSolver->node[iPoint]->GetSensitivity(2);
          }
          SurfAdj_file << endl;
        }
    }
  }
  
  SurfAdj_file.close();
  
#else
  int rank, iProcessor, nProcessor;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  unsigned short nDim = geometry->GetnDim(), iMarker;
  su2double *Solution, *Coord;
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
  
  SU2_MPI::Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  su2double *Buffer_Send_Coord_x = new su2double[MaxLocalVertex_Surface];
  su2double *Buffer_Send_Coord_y= new su2double[MaxLocalVertex_Surface];
  su2double *Buffer_Send_Coord_z= new su2double[MaxLocalVertex_Surface];
  unsigned long *Buffer_Send_GlobalPoint= new unsigned long[MaxLocalVertex_Surface];
  su2double *Buffer_Send_Sensitivity= new su2double[MaxLocalVertex_Surface];
  su2double *Buffer_Send_PsiRho= new su2double[MaxLocalVertex_Surface];
  su2double *Buffer_Send_Phi_x= new su2double[MaxLocalVertex_Surface];
  su2double *Buffer_Send_Phi_y= new su2double[MaxLocalVertex_Surface];
  su2double *Buffer_Send_Phi_z= new su2double[MaxLocalVertex_Surface];
  su2double *Buffer_Send_PsiE= new su2double[MaxLocalVertex_Surface];

  su2double *Buffer_Send_Sens_x = NULL, *Buffer_Send_Sens_y = NULL, *Buffer_Send_Sens_z = NULL;

  if (config->GetDiscrete_Adjoint()){
    Buffer_Send_Sens_x = new su2double[MaxLocalVertex_Surface];
    Buffer_Send_Sens_y = new su2double[MaxLocalVertex_Surface];
    if (nDim == 3){
      Buffer_Send_Sens_z = new su2double[MaxLocalVertex_Surface];
    }
  }
  
  nVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          Solution = AdjSolver->node[iPoint]->GetSolution();
          //Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->node[iPoint]->GetCoord();
          //d = AdjSolver->node[iPoint]->GetForceProj_Vector();
          Buffer_Send_GlobalPoint[nVertex_Surface] = geometry->node[iPoint]->GetGlobalIndex();
          Buffer_Send_Coord_x[nVertex_Surface] = Coord[0];
          Buffer_Send_Coord_y[nVertex_Surface] = Coord[1];
          Buffer_Send_Sensitivity[nVertex_Surface] =  AdjSolver->GetCSensitivity(iMarker, iVertex);
          Buffer_Send_PsiRho[nVertex_Surface] = Solution[0];
          Buffer_Send_Phi_x[nVertex_Surface] = Solution[1];
          Buffer_Send_Phi_y[nVertex_Surface] = Solution[2];
          if (nDim == 2) Buffer_Send_PsiE[nVertex_Surface] = Solution[3];
          if (nDim == 3) {
            Buffer_Send_Coord_z[nVertex_Surface] = Coord[2];
            Buffer_Send_Phi_z[nVertex_Surface] = Solution[3];
            Buffer_Send_PsiE[nVertex_Surface] = Solution[4];
          }
          if (config->GetDiscrete_Adjoint()){
            Buffer_Send_Sens_x[nVertex_Surface] = AdjSolver->node[iPoint]->GetSensitivity(0);
            Buffer_Send_Sens_y[nVertex_Surface] = AdjSolver->node[iPoint]->GetSensitivity(1);
            if (nDim == 3){
              Buffer_Send_Sens_z[nVertex_Surface] = AdjSolver->node[iPoint]->GetSensitivity(2);
            }
          }
          
          /*--- If US system, the output should be in inches ---*/
          
          if (config->GetSystemMeasurements() == US) {
            Buffer_Send_Coord_x[nVertex_Surface] *= 12.0;
            Buffer_Send_Coord_y[nVertex_Surface] *= 12.0;
            if (nDim == 3) Buffer_Send_Coord_z[nVertex_Surface] *= 12.0;
          }
          
          nVertex_Surface++;
        }
      }
  
  su2double *Buffer_Receive_Coord_x = NULL, *Buffer_Receive_Coord_y = NULL, *Buffer_Receive_Coord_z = NULL, *Buffer_Receive_Sensitivity = NULL,
  *Buffer_Receive_PsiRho = NULL, *Buffer_Receive_Phi_x = NULL, *Buffer_Receive_Phi_y = NULL, *Buffer_Receive_Phi_z = NULL,
  *Buffer_Receive_PsiE = NULL, *Buffer_Receive_Sens_x = NULL, *Buffer_Receive_Sens_y = NULL, *Buffer_Receive_Sens_z = NULL;
  unsigned long *Buffer_Receive_GlobalPoint = NULL;
  
  if (rank == MASTER_NODE) {
    Buffer_Receive_Coord_x = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Coord_y = new su2double [nProcessor*MaxLocalVertex_Surface];
    if (nDim == 3) Buffer_Receive_Coord_z = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Sensitivity = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_PsiRho = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Phi_x = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Phi_y = new su2double [nProcessor*MaxLocalVertex_Surface];
    if (nDim == 3) Buffer_Receive_Phi_z = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_PsiE = new su2double [nProcessor*MaxLocalVertex_Surface];
    if (config->GetDiscrete_Adjoint()){
      Buffer_Receive_Sens_x = new su2double[nProcessor*MaxLocalVertex_Surface];
      Buffer_Receive_Sens_y = new su2double[nProcessor*MaxLocalVertex_Surface];
      if (nDim == 3){
        Buffer_Receive_Sens_z = new su2double[nProcessor*MaxLocalVertex_Surface];
      }
    }
  }
  
  nBuffer_Scalar = MaxLocalVertex_Surface;
  
  /*--- Send the information to the Master node ---*/
  SU2_MPI::Gather(Buffer_Send_Coord_x, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Coord_x, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Coord_y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Coord_y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) SU2_MPI::Gather(Buffer_Send_Coord_z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Coord_z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_GlobalPoint, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalPoint, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Sensitivity, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Sensitivity, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_PsiRho, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_PsiRho, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Phi_x, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Phi_x, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Phi_y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Phi_y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) SU2_MPI::Gather(Buffer_Send_Phi_z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Phi_z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_PsiE, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_PsiE, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (config->GetDiscrete_Adjoint()){
    SU2_MPI::Gather(Buffer_Send_Sens_x, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Sens_x, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_Sens_y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Sens_y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    if (nDim == 3){
      SU2_MPI::Gather(Buffer_Send_Sens_z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Sens_z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    }
  }
  
  /*--- The master node is the one who writes the surface files ---*/
  if (rank == MASTER_NODE) {
    unsigned long iVertex, GlobalPoint, position;
    char cstr[200], buffer[50];
    ofstream SurfAdj_file;
    string filename = config->GetSurfAdjCoeff_FileName();
        
    /*--- Write file name with extension if unsteady ---*/
    strcpy (cstr, filename.c_str());
    
    if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
      if (SU2_TYPE::Int(val_iZone) < 10) SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(val_iZone));
      if ((SU2_TYPE::Int(val_iZone) >= 10) && (SU2_TYPE::Int(val_iZone) < 100)) SPRINTF (buffer, "_000%d.csv", SU2_TYPE::Int(val_iZone));
      if ((SU2_TYPE::Int(val_iZone) >= 100) && (SU2_TYPE::Int(val_iZone) < 1000)) SPRINTF (buffer, "_00%d.csv", SU2_TYPE::Int(val_iZone));
      if ((SU2_TYPE::Int(val_iZone) >= 1000) && (SU2_TYPE::Int(val_iZone) < 10000)) SPRINTF (buffer, "_0%d.csv", SU2_TYPE::Int(val_iZone));
      if (SU2_TYPE::Int(val_iZone) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(val_iZone));
      
    } else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
      if ((SU2_TYPE::Int(iExtIter) >= 0) && (SU2_TYPE::Int(iExtIter) < 10)) SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 10) && (SU2_TYPE::Int(iExtIter) < 100)) SPRINTF (buffer, "_000%d.csv", SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 100) && (SU2_TYPE::Int(iExtIter) < 1000)) SPRINTF (buffer, "_00%d.csv", SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.csv", SU2_TYPE::Int(iExtIter));
      if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(iExtIter));
    }
    else
      SPRINTF (buffer, ".csv");
    
    strcat (cstr, buffer);
    SurfAdj_file.open(cstr, ios::out);
    SurfAdj_file.precision(15);
    
    /*--- Write the 2D surface flow coefficient file ---*/
    if (geometry->GetnDim() == 2) {
      
      SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"";
      if (config->GetDiscrete_Adjoint()){
        SurfAdj_file << ",\" x_Sens\",\"y_Sens\"";
      }
      SurfAdj_file << endl;

      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
        for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
          
          position = iProcessor*MaxLocalVertex_Surface+iVertex;
          GlobalPoint = Buffer_Receive_GlobalPoint[position];
          
          SurfAdj_file << scientific << GlobalPoint <<
          ", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
          ", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] <<
          ", " << Buffer_Receive_PsiE[position] << ", " << Buffer_Receive_Coord_x[position] <<
          ", "<< Buffer_Receive_Coord_y[position];
          if (config->GetDiscrete_Adjoint()){
            SurfAdj_file << ", " << Buffer_Receive_Sens_x[position] << ", " << Buffer_Receive_Sens_y[position];
          }
          SurfAdj_file << endl;
        }
    }
    
    /*--- Write the 3D surface flow coefficient file ---*/
    if (geometry->GetnDim() == 3) {
      
      SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"" << endl;
      if (config->GetDiscrete_Adjoint()){
        SurfAdj_file << ",\"x_Sens\",\"y_Sens\",\"z_Sens\"";
      }
      SurfAdj_file << endl;

      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
        for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
          position = iProcessor*MaxLocalVertex_Surface+iVertex;
          GlobalPoint = Buffer_Receive_GlobalPoint[position];
          
          SurfAdj_file << scientific << GlobalPoint <<
          ", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
          ", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] << ", " << Buffer_Receive_Phi_z[position] <<
          ", " << Buffer_Receive_PsiE[position] <<", "<< Buffer_Receive_Coord_x[position] <<
          ", "<< Buffer_Receive_Coord_y[position] <<", "<< Buffer_Receive_Coord_z[position];
          if (config->GetDiscrete_Adjoint()){
            SurfAdj_file << ", " << Buffer_Receive_Sens_x[position] << ", " << Buffer_Receive_Sens_y[position] << ", " << Buffer_Receive_Sens_z[position];
          }
          SurfAdj_file << endl;
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
    if (config->GetDiscrete_Adjoint()){
      delete [] Buffer_Receive_Sens_x;
      delete [] Buffer_Receive_Sens_y;
      if (nDim == 3){
        delete [] Buffer_Receive_Sens_z;
      }
    }
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
  if (Buffer_Send_Sens_x != NULL) delete [] Buffer_Send_Sens_x;
  if (Buffer_Send_Sens_y != NULL) delete [] Buffer_Send_Sens_y;
  if (Buffer_Send_Sens_z != NULL) delete [] Buffer_Send_Sens_z;
  
  SurfAdj_file.close();
  
#endif
}

void COutput::MergeConnectivity(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  /*--- Flags identifying the types of files to be written. ---*/
  
  bool Wrt_Vol = config->GetWrt_Vol_Sol();
  bool Wrt_Srf = config->GetWrt_Srf_Sol();
  
  /*--- Merge connectivity for each type of element (excluding halos). Note
   that we only need to merge the connectivity once, as it does not change
   during computation. Check whether the base file has been written. ---*/
  
    /*--- Merge volumetric grid. ---*/
    
    if (Wrt_Vol) {
      
      if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Tria != 0))
        cout <<"Merging volumetric triangle grid connectivity." << endl;
      MergeVolumetricConnectivity(config, geometry, TRIANGLE    );
      
      if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Quad != 0))
        cout <<"Merging volumetric quadrilateral grid connectivity." << endl;
      MergeVolumetricConnectivity(config, geometry, QUADRILATERAL   );
      
      if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Tetr != 0))
        cout <<"Merging volumetric tetrahedron grid connectivity." << endl;
      MergeVolumetricConnectivity(config, geometry, TETRAHEDRON );
      
      if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Hexa != 0))
        cout <<"Merging volumetric hexahedron grid connectivity." << endl;
      MergeVolumetricConnectivity(config, geometry, HEXAHEDRON  );
      
      if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Pris != 0))
        cout <<"Merging volumetric prism grid connectivity." << endl;
      MergeVolumetricConnectivity(config, geometry, PRISM       );
      
      if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Pyra != 0))
        cout <<"Merging volumetric pyramid grid connectivity." << endl;
      MergeVolumetricConnectivity(config, geometry, PYRAMID     );
      
    }
    
    /*--- Merge surface grid. ---*/
    
    if (Wrt_Srf) {
      
      if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Line != 0))
        cout <<"Merging surface line grid connectivity." << endl;
      MergeSurfaceConnectivity(config, geometry, LINE);
      
      if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_BoundTria != 0))
        cout <<"Merging surface triangle grid connectivity." << endl;
      MergeSurfaceConnectivity(config, geometry, TRIANGLE);
      
      if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_BoundQuad != 0))
        cout <<"Merging surface quadrilateral grid connectivity." << endl;
      MergeSurfaceConnectivity(config, geometry, QUADRILATERAL);
      
    }
    
    /*--- Update total number of volume elements after merge. ---*/
    
    nGlobal_Elem = nGlobal_Tria + nGlobal_Quad + nGlobal_Tetr +
    nGlobal_Hexa + nGlobal_Pyra + nGlobal_Pris;
    
    /*--- Update total number of surface elements after merge. ---*/
    
    nSurf_Elem = nGlobal_Line + nGlobal_BoundTria + nGlobal_BoundQuad;

}

void COutput::MergeCoordinates(CConfig *config, CGeometry *geometry) {
  
  /*--- Local variables needed on all processors ---*/
  
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint;
  
#ifndef HAVE_MPI
  
  /*--- In serial, the single process has access to all geometry, so simply
   load the coordinates into the data structure. ---*/
  
  unsigned short iMarker;
  unsigned long iVertex, nTotalPoints = 0;
  int SendRecv;
  
  /*--- First, create a structure to locate any periodic halo nodes ---*/
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
            (SendRecv < 0)) {
          Local_Halo[iPoint] = false;
        }
      }
      
    }
  }
  
  /*--- Total number of points in the mesh (this might include periodic points). ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (!Local_Halo[iPoint]) nTotalPoints++;
  
  nGlobal_Poin = nTotalPoints;
  nGlobal_Doma = geometry->GetnPointDomain();
  
  /*--- Allocate the coordinates data structure. ---*/
  
  Coords = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Coords[iDim] = new su2double[nGlobal_Poin];
  }
  
  /*--- Loop over the mesh to collect the coords of the local points ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node). 
     Sort by the global index, even in serial there is a renumbering (e.g. RCM). ---*/
    
    if (!Local_Halo[iPoint]) {
      
      /*--- Retrieve the current coordinates at this node. ---*/
      
      unsigned long iGlobal_Index = geometry->node[iPoint]->GetGlobalIndex();
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Coords[iDim][iGlobal_Index] = geometry->node[iPoint]->GetCoord(iDim);
        
        /*--- If US system, the output should be in inches ---*/
        
        if ((config->GetSystemMeasurements() == US) && (config->GetKind_SU2() != SU2_DEF)) {
          Coords[iDim][iGlobal_Index] *= 12.0;
        }
        
      }
      
    }
  }

  
  delete [] Local_Halo;
  
#else
  
  /*--- MPI preprocessing ---*/
  int iProcessor, nProcessor, rank;
  unsigned long jPoint;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;
  
  /*--- Local variables needed for merging the geometry with MPI. ---*/
  
  unsigned long iVertex, iMarker;
  unsigned long Buffer_Send_nPoin[1], *Buffer_Recv_nPoin = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  
  if (rank == MASTER_NODE) Buffer_Recv_nPoin = new unsigned long[nProcessor];
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/
  
  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        
        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }
    
    /*--- Sum total number of nodes that belong to the domain ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;
  }
  Buffer_Send_nPoin[0] = nLocalPoint;
  
  /*--- Communicate the total number of nodes on this domain. ---*/
  
  SU2_MPI::Gather(&Buffer_Send_nPoin, 1, MPI_UNSIGNED_LONG,
             Buffer_Recv_nPoin, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  
  if (rank == MASTER_NODE) {
    nGlobal_Doma = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Doma += Buffer_Recv_nPoin[iProcessor];
    }
  }
  nBuffer_Scalar = MaxLocalPoint;
  
  /*--- Send and Recv buffers. ---*/
  
  su2double *Buffer_Send_X = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_X = NULL;
  
  su2double *Buffer_Send_Y = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Y = NULL;
  
  su2double *Buffer_Send_Z = NULL, *Buffer_Recv_Z = NULL;
  if (nDim == 3) Buffer_Send_Z = new su2double[MaxLocalPoint];
  
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Prepare the receive buffers in the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    
    Buffer_Recv_X = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Y = new su2double[nProcessor*MaxLocalPoint];
    if (nDim == 3) Buffer_Recv_Z = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];
    
    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoin[iProcessor];
    }
    Coords = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Coords[iDim] = new su2double[nGlobal_Poin];
    }
  }
  
  /*--- Main communication routine. Loop over each coordinate and perform
   the MPI comm. Temporary 1-D buffers are used to send the coordinates at
   all nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/
  
  /*--- Loop over this partition to collect the coords of the local points. ---*/
  su2double *Coords_Local; jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check for halos and write only if requested ---*/
    if (!Local_Halo[iPoint] || Wrt_Halo) {
      
      /*--- Retrieve local coordinates at this node. ---*/
      Coords_Local = geometry->node[iPoint]->GetCoord();
      
      /*--- Load local coords into the temporary send buffer. ---*/
      Buffer_Send_X[jPoint] = Coords_Local[0];
      Buffer_Send_Y[jPoint] = Coords_Local[1];
      if (nDim == 3) Buffer_Send_Z[jPoint] = Coords_Local[2];
      
      /*--- If US system, the output should be in inches ---*/
      
      if ((config->GetSystemMeasurements() == US) && (config->GetKind_SU2() != SU2_DEF)) {
        Buffer_Send_X[jPoint] *= 12.0;
        Buffer_Send_Y[jPoint] *= 12.0;
        if (nDim == 3) Buffer_Send_Z[jPoint] *= 12.0;
      }
      
      /*--- Store the global index for this local node. ---*/
      Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
      
      /*--- Increment jPoint as the counter. We need this because iPoint
       may include halo nodes that we skip over during this loop. ---*/
      jPoint++;
    }
  }

  /*--- Gather the coordinate data on the master node using MPI. ---*/

  SU2_MPI::Gather(Buffer_Send_X, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_X, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) {
    SU2_MPI::Gather(Buffer_Send_Z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  }
  SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

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
  
  delete [] Local_Halo;
  delete [] Buffer_Send_X;
  delete [] Buffer_Send_Y;
  if (Buffer_Send_Z != NULL) delete [] Buffer_Send_Z;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_X;
    delete [] Buffer_Recv_Y;
    if (Buffer_Recv_Z != NULL)  delete [] Buffer_Recv_Z;
    delete [] Buffer_Recv_GlobalIndex;
    delete [] Buffer_Recv_nPoin;
  }
  
#endif
  
}

void COutput::MergeVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  int iProcessor;
  unsigned short NODES_PER_ELEMENT;
  unsigned long iPoint, iNode, jNode;
  unsigned long iElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;
  
  unsigned long iVertex, iMarker;
  unsigned long jElem;
  int SendRecv, RecvFrom;
  
  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0;
  unsigned long MaxLocalElem = 0, iGlobal_Index, jPoint, kPoint;
  
  bool Wrt_Halo = config->GetWrt_Halo();
  bool *Write_Elem = NULL, notPeriodic, notHalo, addedPeriodic;

  int *Conn_Elem = NULL;

  int rank = MASTER_NODE;
  int size = SINGLE_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  
  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/
  
  switch (Elem_Type) {
    case TRIANGLE:
      nLocalElem = geometry->GetnElemTria();
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
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
    case PRISM:
      nLocalElem = geometry->GetnElemPris();
      NODES_PER_ELEMENT = N_POINTS_PRISM;
      break;
    case PYRAMID:
      nLocalElem = geometry->GetnElemPyra();
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      cout << "Error: Unrecognized element type \n";
      exit(EXIT_FAILURE); break;
  }
  
  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/
  
  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[size];
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalElem, &MaxLocalElem, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nElem, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nElem, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalElem = nLocalElem;
  Buffer_Recv_nElem[0] = Buffer_Send_nElem[0];
#endif
  
  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;
  
  /*--- Send and Recv buffers ---*/
  
  unsigned long *Buffer_Send_Elem = new unsigned long[nBuffer_Scalar];
  unsigned long *Buffer_Recv_Elem = NULL;
  
  unsigned short *Buffer_Send_Halo = new unsigned short[MaxLocalElem];
  unsigned short *Buffer_Recv_Halo = NULL;
  
  /*--- Prepare the receive buffers on the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new unsigned long[size*nBuffer_Scalar];
    Buffer_Recv_Halo = new unsigned short[size*MaxLocalElem];
    Conn_Elem = new int[size*MaxLocalElem*NODES_PER_ELEMENT];
  }
  
  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that we part of the original domain (we want these to be in the
   output files). ---*/
  
  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
            (SendRecv < 0)) {
          Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
        }
      }
    }
  }
  
  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/
  
  /*--- First communicate the number of points that each rank has found ---*/
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];
  
  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif
  
  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];
  
  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }
  
  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++) Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iGlobal_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- We need to keep one copy of overlapping halo cells. ---*/
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));
        
        /*--- We want to keep the periodic nodes that were part of the original domain ---*/
        notPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                       (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
                       (SendRecv < 0));
        
        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (iGlobal_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          kPoint = (iProcessor+1)*maxAddedPeriodic;
        }
        
        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }
      }
    }
  }
  
  /*--- Loop over all elements in this partition and load the
   elements of the current type into the buffer to be sent to
   the master node. ---*/
  
  jNode = 0; jElem = 0;
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if (geometry->elem[iElem]->GetVTK_Type() == Elem_Type) {
      
      /*--- Loop over all nodes in this element and load the
       connectivity into the send buffer. ---*/
      
      Buffer_Send_Halo[jElem] = false;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        
        /*--- Store the global index values directly. ---*/
        
        iPoint = geometry->elem[iElem]->GetNode(iNode);
        Buffer_Send_Elem[jNode] = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- Check if this is a halo node. If so, flag this element
         as a halo cell. We will use this later to sort and remove
         any duplicates from the connectivity list. ---*/
        
        if (Local_Halo[iPoint]) {
          Buffer_Send_Halo[jElem] = true;
        }
        
        /*--- Increment jNode as the counter. We need this because iElem
         may include other elements that we skip over during this loop. ---*/
        
        jNode++;
      }
      jElem++;
    }
  }
  
  /*--- Gather the element connectivity information. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, Buffer_Recv_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Elem[iPoint] = Buffer_Send_Elem[iPoint];
  for (iPoint = 0; iPoint < MaxLocalElem; iPoint++) Buffer_Recv_Halo[iPoint] = Buffer_Send_Halo[iPoint];
#endif
  
  /*--- The master node unpacks and sorts the connectivity. ---*/
  
  if (rank == MASTER_NODE) {
    
    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/
    
    Write_Elem = new bool[size*MaxLocalElem];
    for (iElem = 0; iElem < size*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }
    
    /*--- Remove the rind layer from the solution only if requested ---*/
    
    if (!Wrt_Halo) {
      
      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/
      
      kElem = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
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
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
        
        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {
          
          /*--- Increment total count for this element type ---*/
          nElem_Total++;
          
          /*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/
          
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = (int)Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
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
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] Local_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }
  
  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/
  
  if (rank == MASTER_NODE) {
    switch (Elem_Type) {
      case TRIANGLE:
        nGlobal_Tria = nElem_Total;
        if (nGlobal_Tria > 0) Conn_Tria = Conn_Elem;
        break;
      case QUADRILATERAL:
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
      case PRISM:
        nGlobal_Pris = nElem_Total;
        if (nGlobal_Pris > 0) Conn_Pris = Conn_Elem;
        break;
      case PYRAMID:
        nGlobal_Pyra = nElem_Total;
        if (nGlobal_Pyra > 0) Conn_Pyra = Conn_Elem;
        break;
      default:
        cout << "Error: Unrecognized element type \n";
        exit(EXIT_FAILURE); break;
    }
  }
  
}

void COutput::MergeSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  unsigned short NODES_PER_ELEMENT;
  
  unsigned short iMarker;
  unsigned long iPoint, iNode, jNode;
  unsigned long iElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;
  
  int iProcessor;
  unsigned long jElem;
  
  unsigned long iVertex;
  
  int SendRecv, RecvFrom;
  
  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0;
  unsigned long MaxLocalElem = 0, iGlobal_Index, jPoint, kPoint;
  
  bool Wrt_Halo = config->GetWrt_Halo();
  bool *Write_Elem = NULL, notPeriodic, notHalo, addedPeriodic;

  
  int *Conn_Elem = NULL;

  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
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
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      cout << "Error: Unrecognized element type \n";
      exit(EXIT_FAILURE); break;
  }
  
  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/
  
  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[size];
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalElem, &MaxLocalElem, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nElem, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nElem, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalElem = nLocalElem;
  Buffer_Recv_nElem[0] = Buffer_Send_nElem[0];
#endif
  
  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;
  
  /*--- Send and Recv buffers ---*/
  
  unsigned long *Buffer_Send_Elem = new unsigned long[nBuffer_Scalar];
  unsigned long *Buffer_Recv_Elem = NULL;
  
  unsigned short *Buffer_Send_Halo = new unsigned short[MaxLocalElem];
  unsigned short *Buffer_Recv_Halo = NULL;
  
  /*--- Prepare the receive buffers on the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new unsigned long[size*nBuffer_Scalar];
    Buffer_Recv_Halo = new unsigned short[size*MaxLocalElem];
    Conn_Elem = new int[size*MaxLocalElem*NODES_PER_ELEMENT];
  }
  
  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that we part of the original domain (we want these to be in the
   output files). ---*/
  
  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
            (SendRecv < 0)) {
          Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
        }
      }
    }
  }
  
  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/
  
  /*--- First communicate the number of points that each rank has found ---*/
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];
  
  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif
  
  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];
  
  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }
  
  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++) Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iGlobal_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- We need to keep one copy of overlapping halo cells. ---*/
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));
        
        /*--- We want to keep the periodic nodes that were part of the original domain ---*/
        notPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                       (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
                       (SendRecv < 0));
        
        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (iGlobal_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          kPoint = (iProcessor+1)*maxAddedPeriodic;
        }
        
        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
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
      for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
        
        if (geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
          
          /*--- Loop over all nodes in this element and load the
           connectivity into the send buffer. ---*/
          
          Buffer_Send_Halo[jElem] = false;
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            
            /*--- Store the global index values directly. ---*/
            
            iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
            Buffer_Send_Elem[jNode] = geometry->node[iPoint]->GetGlobalIndex();
            
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

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, Buffer_Recv_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Elem[iPoint] = Buffer_Send_Elem[iPoint];
  for (iPoint = 0; iPoint < MaxLocalElem; iPoint++) Buffer_Recv_Halo[iPoint] = Buffer_Send_Halo[iPoint];
#endif
  
  /*--- The master node unpacks and sorts the connectivity. ---*/
  
  if (rank == MASTER_NODE) {
    
    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/
    
    Write_Elem = new bool[size*MaxLocalElem];
    for (iElem = 0; iElem < size*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }
    
    /*--- Remove the rind layer from the solution only if requested ---*/
    
    if (!Wrt_Halo) {
      
      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/
      
      kElem = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
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
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
        
        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {
          
          /*--- Increment total count for this element type ---*/
          nElem_Total++;
          
          /*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/
          
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = (int)Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
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
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] Local_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }
  
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
      case QUADRILATERAL:
        nGlobal_BoundQuad = nElem_Total;
        if (nGlobal_BoundQuad > 0) Conn_BoundQuad = Conn_Elem;
        break;
      default:
        cout << "Error: Unrecognized element type \n";
        exit(EXIT_FAILURE); break;
    }
  }
  
}

void COutput::MergeSolution(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar = 0, jVar = 0, FirstIndex = NONE, SecondIndex = NONE, ThirdIndex = NONE;
  unsigned short nVar_First = 0, nVar_Second = 0, nVar_Third = 0;
  unsigned short iVar_GridVel = 0, iVar_PressCp = 0, iVar_Density = 0, iVar_Lam = 0, iVar_MachMean = 0,
  iVar_ViscCoeffs = 0, iVar_Sens = 0, iVar_Extra = 0, iVar_Eddy = 0, iVar_Sharp = 0,
  iVar_FEA_Stress = 0, iVar_FEA_Stress_3D = 0, iVar_FEA_Extra = 0, iVar_SensDim = 0;
  
  unsigned long iPoint = 0, jPoint = 0, iVertex = 0, iMarker = 0;
  su2double Gas_Constant, Mach2Vel, Mach_Motion, RefDensity, RefPressure = 0.0, factor = 0.0;
  
  su2double *Aux_Frict = NULL, *Aux_Heat = NULL, *Aux_yPlus = NULL, *Aux_Sens = NULL;
  
  unsigned short CurrentIndex;
  int *Local_Halo;
  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;

  int iProcessor;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
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
  unsigned short nDim = geometry->GetnDim();
  su2double RefAreaCoeff = config->GetRefAreaCoeff();
  su2double Gamma = config->GetGamma();
  su2double RefVel2, *Normal, Area;
  
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
    case EULER : case NAVIER_STOKES: FirstIndex = FLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case RANS : FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; if (transition) ThirdIndex=TRANS_SOL; else ThirdIndex = NONE; break;
    case POISSON_EQUATION: FirstIndex = POISSON_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case WAVE_EQUATION: FirstIndex = WAVE_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case HEAT_EQUATION: FirstIndex = HEAT_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case LINEAR_ELASTICITY: FirstIndex = FEA_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case ADJ_EULER : case ADJ_NAVIER_STOKES : FirstIndex = ADJFLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case ADJ_RANS : FirstIndex = ADJFLOW_SOL; if (config->GetFrozen_Visc()) SecondIndex = NONE; else SecondIndex = ADJTURB_SOL; ThirdIndex = NONE; break;
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: FirstIndex = ADJFLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case DISC_ADJ_RANS: FirstIndex = ADJFLOW_SOL; SecondIndex = ADJTURB_SOL; ThirdIndex = NONE; break;
    default: SecondIndex = NONE; ThirdIndex = NONE; break;
  }
  
  nVar_First = solver[FirstIndex]->GetnVar();
  if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
  if (ThirdIndex != NONE) nVar_Third = solver[ThirdIndex]->GetnVar();
  nVar_Consv = nVar_First + nVar_Second + nVar_Third;
  nVar_Total = nVar_Consv;
  
  if (!config->GetLow_MemoryOutput()) {
    
    /*--- Add the limiters ---*/
    
    if (config->GetWrt_Limiters()) nVar_Total += nVar_Consv;
    
    /*--- Add the residuals ---*/
    
    if (config->GetWrt_Residuals()) nVar_Total += nVar_Consv;
    
    /*--- Add the grid velocity to the restart file for the unsteady adjoint ---*/
    
    if (grid_movement) {
      iVar_GridVel = nVar_Total;
      if (geometry->GetnDim() == 2) nVar_Total += 2;
      else if (geometry->GetnDim() == 3) nVar_Total += 3;
    }
    
    /*--- Add density to the restart file ---*/

    if ((config->GetKind_Regime() == FREESURFACE)) {
      iVar_Density = nVar_Total; nVar_Total += 1;
    }
    
    /*--- Add Pressure, Temperature, Cp, Mach to the restart file ---*/

    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      iVar_PressCp = nVar_Total; nVar_Total += 3;
      iVar_MachMean = nVar_Total; nVar_Total += 1;
    }
    
    /*--- Add Laminar Viscosity, Skin Friction, Heat Flux, & yPlus to the restart file ---*/

    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      iVar_Lam = nVar_Total; nVar_Total += 1;
      iVar_ViscCoeffs = nVar_Total; nVar_Total += 3;
    }
    
    /*--- Add Eddy Viscosity to the restart file ---*/

    if (Kind_Solver == RANS) {
      iVar_Eddy = nVar_Total; nVar_Total += 1;
    }
    
    /*--- Add Sharp edges to the restart file ---*/

    if (config->GetWrt_SharpEdges()) {
      if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
        iVar_Sharp = nVar_Total; nVar_Total += 1;
      }
    }
    
    //if (Kind_Solver == POISSON_EQUATION) {
    //  iVar_EF = nVar_Total; nVar_Total += geometry->GetnDim();
    //}
    
    if (( Kind_Solver == ADJ_EULER              ) || ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
        ( Kind_Solver == ADJ_RANS               )) {
      iVar_Sens   = nVar_Total; nVar_Total += 2;
    }
    
    if (Kind_Solver == LINEAR_ELASTICITY) {
      iVar_FEA_Stress  = nVar_Total; nVar_Total += 3;
	  if (geometry->GetnDim() == 3) {iVar_FEA_Stress_3D = nVar_Total; nVar_Total += 3;}
      iVar_FEA_Extra = nVar_Total; nVar_Total += 2;
    }

    if ((Kind_Solver == DISC_ADJ_EULER)         ||
        (Kind_Solver == DISC_ADJ_NAVIER_STOKES) ||
         (Kind_Solver == DISC_ADJ_RANS)){
      iVar_Sens    = nVar_Total; nVar_Total += 1;
      iVar_SensDim = nVar_Total; nVar_Total += nDim;
    }
    
    if (config->GetExtraOutput()) {
      if (Kind_Solver == RANS) {
        iVar_Extra  = nVar_Total; nVar_Extra  = solver[TURB_SOL]->GetnOutputVariables(); nVar_Total += nVar_Extra;
      }
    }
    
  }
  
  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/
  
  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        
        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }
    
    /*--- Sum total number of nodes that belong to the domain ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;
    
  }
  Buffer_Send_nPoint[0] = nLocalPoint;
  
  /*--- Each processor sends its local number of nodes to the master. ---*/
  
  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[size];
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalPoint = nLocalPoint;
  Buffer_Recv_nPoint[0] = Buffer_Send_nPoint[0];
#endif
  
  nBuffer_Scalar = MaxLocalPoint;
  
  /*--- Send and Recv buffers. ---*/
  
  su2double *Buffer_Send_Var = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Var = NULL;
  
  su2double *Buffer_Send_Res = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Res = NULL;
  
  su2double *Buffer_Send_Vol = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Vol = NULL;
  
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Auxiliary vectors for surface coefficients ---*/
  
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    Aux_Frict = new su2double[geometry->GetnPoint()];
    Aux_Heat  = new su2double[geometry->GetnPoint()];
    Aux_yPlus = new su2double[geometry->GetnPoint()];
  }
  
  if ((Kind_Solver == ADJ_EULER) ||
      (Kind_Solver == ADJ_NAVIER_STOKES) ||
      (Kind_Solver == ADJ_RANS)  ||
      (Kind_Solver == DISC_ADJ_EULER) ||
      (Kind_Solver == DISC_ADJ_NAVIER_STOKES) ||
      (Kind_Solver == DISC_ADJ_RANS)) {
    Aux_Sens = new su2double[geometry->GetnPoint()];
  }
  
  /*--- Prepare the receive buffers in the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    
    Buffer_Recv_Var = new su2double[size*MaxLocalPoint];
    Buffer_Recv_Res = new su2double[size*MaxLocalPoint];
    Buffer_Recv_Vol = new su2double[size*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[size*MaxLocalPoint];
    
    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
    }
    Data = new su2double*[nVar_Total];
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      Data[iVar] = new su2double[nGlobal_Poin];
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
      
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Get this variable into the temporary send buffer. ---*/
        
        Buffer_Send_Var[jPoint] = solver[CurrentIndex]->node[iPoint]->GetSolution(jVar);
        
        if (!config->GetLow_MemoryOutput()) {
          
          if (config->GetWrt_Limiters()) {
            Buffer_Send_Vol[jPoint] = solver[CurrentIndex]->node[iPoint]->GetLimiter_Primitive(jVar);
          }
          
          if (config->GetWrt_Residuals()) {
            Buffer_Send_Res[jPoint] = solver[CurrentIndex]->LinSysRes.GetBlock(iPoint, jVar);
          }
          
        }
        
        /*--- Only send/recv the volumes & global indices during the first loop ---*/
        
        if (iVar == 0) {
          Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
        }
        
        jPoint++;
        
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    
#ifdef HAVE_MPI
    SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif
    if (!config->GetLow_MemoryOutput()) {
      
      if (config->GetWrt_Limiters()) {
#ifdef HAVE_MPI
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif
      }
      
      if (config->GetWrt_Residuals()) {
#ifdef HAVE_MPI
        SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
#endif
      }
      
    }
    
    if (iVar == 0) {
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_GlobalIndex[iPoint] = Buffer_Send_GlobalIndex[iPoint];
#endif
    }
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    
    if (rank == MASTER_NODE) {
      jPoint = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          
          if (!config->GetLow_MemoryOutput()) {
            
            if (config->GetWrt_Limiters()) {
              Data[iVar+nVar_Consv][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            }
            
            if (config->GetWrt_Residuals()) {
              unsigned short ExtraIndex;
              ExtraIndex = nVar_Consv;
              if (config->GetWrt_Limiters()) ExtraIndex = 2*nVar_Consv;
              Data[iVar+ExtraIndex][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            }
            
          }
          
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
    
  }
  
  if (!config->GetLow_MemoryOutput()) {
    
    /*--- Additional communication routine for the grid velocity. Note that
     we are reusing the same temporary buffers from above for efficiency.
     Also, in the future more routines like this could be used to write
     an arbitrary number of additional variables to the file. ---*/
    
    if (grid_movement) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0; su2double *Grid_Vel;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the three grid velocity components. ---*/
          
          Grid_Vel = geometry->node[iPoint]->GetGridVel();
          Buffer_Send_Var[jPoint] = Grid_Vel[0];
          Buffer_Send_Res[jPoint] = Grid_Vel[1];
          if (geometry->GetnDim() == 3) Buffer_Send_Vol[jPoint] = Grid_Vel[2];
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (geometry->GetnDim() == 3) {
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      }
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      if (geometry->GetnDim() == 3) {
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
      }
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_GridVel;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
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
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the pressure and mach variables. ---*/
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc();
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_Density;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
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
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the pressure, Cp, and mach variables. ---*/
          
          if (compressible) {
            Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure();
            Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature();
            Buffer_Send_Vol[jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure() - RefPressure)*factor*RefAreaCoeff;
          }
          if (incompressible || freesurface) {
            Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressureInc();
            Buffer_Send_Res[jPoint] = 0.0;
            Buffer_Send_Vol[jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressureInc() - RefPressure)*factor*RefAreaCoeff;
          }
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_PressCp;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
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
    
    /*--- Communicate Mach---*/
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
          
          if (compressible) {
            Buffer_Send_Var[jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/
            solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed();
          }
          if (incompressible || freesurface) {
            Buffer_Send_Var[jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/
            sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref()));
          }
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_MachMean;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
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
    
    /*--- Laminar Viscosity ---*/
    
    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
          
          if (compressible) {
            Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
          }
          if (incompressible || freesurface) {
            Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc();
          }
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_Lam;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    
      /*--- Communicate skin friction, heat transfer, y+ ---*/
      
      /*--- First, loop through the mesh in order to find and store the
       value of the viscous coefficients at any surface nodes. They
       will be placed in an auxiliary vector and then communicated like
       all other volumetric variables. ---*/
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        Aux_Frict[iPoint] = 0.0;
        Aux_Heat[iPoint]  = 0.0;
        Aux_yPlus[iPoint] = 0.0;
      }
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
        if (config->GetMarker_All_Plotting(iMarker) == YES) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            Aux_Frict[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex);
            Aux_Heat[iPoint]  = solver[FLOW_SOL]->GetHeatFlux(iMarker, iVertex);
            Aux_yPlus[iPoint] = solver[FLOW_SOL]->GetYPlus(iMarker, iVertex);
          }
        }
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
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
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_ViscCoeffs;

        for (iProcessor = 0; iProcessor < size; iProcessor++) {
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
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the pressure and mach variables. ---*/
          
          if (compressible) {
            Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosity();
          }
          if (incompressible || freesurface) {
            Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosityInc();
          }
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_Eddy;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
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
    
    if (config->GetWrt_SharpEdges()) {
      
      if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
        
        /*--- Loop over this partition to collect the current variable ---*/
        jPoint = 0;
        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
          
          /*--- Check for halos & write only if requested ---*/
          
          if (!Local_Halo[iPoint] || Wrt_Halo) {
            
            /*--- Load buffers with the pressure and mach variables. ---*/
            
            Buffer_Send_Var[jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance();
            jPoint++;
          }
        }
        
        /*--- Gather the data on the master node. ---*/
        
#ifdef HAVE_MPI
        SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif
        
        /*--- The master node unpacks and sorts this variable by global index ---*/
        
        if (rank == MASTER_NODE) {
          jPoint = 0; iVar = iVar_Sharp;
          for (iProcessor = 0; iProcessor < size; iProcessor++) {
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
    }
    
    /*--- Communicate the surface sensitivity ---*/
    
    if ((Kind_Solver == ADJ_EULER)         ||
        (Kind_Solver == ADJ_NAVIER_STOKES) ||
        (Kind_Solver == ADJ_RANS)          ||
        (Kind_Solver == DISC_ADJ_EULER)    ||
        (Kind_Solver == DISC_ADJ_NAVIER_STOKES) ||
        (Kind_Solver == DISC_ADJ_RANS)) {
      
      /*--- First, loop through the mesh in order to find and store the
       value of the surface sensitivity at any surface nodes. They
       will be placed in an auxiliary vector and then communicated like
       all other volumetric variables. ---*/
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Aux_Sens[iPoint] = 0.0;
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
        if (config->GetMarker_All_Plotting(iMarker) == YES) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
            Area = sqrt (Area);
            Aux_Sens[iPoint] = solver[ADJFLOW_SOL]->GetCSensitivity(iMarker, iVertex)/Area;
          }
        }
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/
          
          Buffer_Send_Var[jPoint] = Aux_Sens[iPoint];
          if ((config->GetKind_ConvNumScheme() == SPACE_CENTERED) && (!config->GetDiscrete_Adjoint()))
            Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetSensor(iPoint);
          if ((config->GetKind_ConvNumScheme() == SPACE_UPWIND) && (!config->GetDiscrete_Adjoint()))
            Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetLimiter(0);
          
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (!config->GetDiscrete_Adjoint())
        SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      if (!config->GetDiscrete_Adjoint())
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_Sens;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar+0][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            if (!config->GetDiscrete_Adjoint())
              Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }

    if ((Kind_Solver == DISC_ADJ_EULER)    ||
        (Kind_Solver == DISC_ADJ_NAVIER_STOKES) ||
        (Kind_Solver == DISC_ADJ_RANS)) {
      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/

          Buffer_Send_Var[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetSensitivity(0);
          Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetSensitivity(1);
          if (nDim == 3)
            Buffer_Send_Vol[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetSensitivity(2);
          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (nDim == 3)
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      if (nDim == 3)
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_SensDim;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar+0][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            if (nDim == 3)
              Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }

    /*--- Communicate the Linear elasticity stresses (2D) ---*/

    if (Kind_Solver == LINEAR_ELASTICITY) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0; su2double **Stress;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the three grid velocity components. ---*/
          
          Stress = solver[FEA_SOL]->node[iPoint]->GetStress();
          /*--- Sigma xx ---*/
          Buffer_Send_Var[jPoint] = Stress[0][0];
          /*--- Sigma yy ---*/
          Buffer_Send_Res[jPoint] = Stress[1][1];
          /*--- Sigma xy ---*/
          Buffer_Send_Vol[jPoint] = Stress[0][1];
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_FEA_Stress;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
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
    
    /*--- Communicate the Linear elasticity stresses (3D) ---*/    
    
    if ((Kind_Solver == LINEAR_ELASTICITY) && (geometry->GetnDim() == 3)) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0; su2double **Stress;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the three grid velocity components. ---*/
          
          Stress = solver[FEA_SOL]->node[iPoint]->GetStress();
          /*--- Sigma zz ---*/
          Buffer_Send_Var[jPoint] = Stress[2][2];
          /*--- Sigma xz ---*/
          Buffer_Send_Res[jPoint] = Stress[0][2];
          /*--- Sigma yz ---*/
          Buffer_Send_Vol[jPoint] = Stress[1][2]; 
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];

#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_FEA_Stress_3D;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
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

    
    /*--- Communicate the Linear elasticity ---*/
    
    if ( Kind_Solver == LINEAR_ELASTICITY ) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
          
          Buffer_Send_Var[jPoint] = solver[FEA_SOL]->node[iPoint]->GetVonMises_Stress();
          Buffer_Send_Res[jPoint] = solver[FEA_SOL]->node[iPoint]->GetFlow_Pressure();
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_FEA_Extra;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
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
          
          if (!Local_Halo[iPoint] || Wrt_Halo) {
            
            /*--- Get this variable into the temporary send buffer. ---*/
            
            if (Kind_Solver == RANS) {
              Buffer_Send_Var[jPoint] = solver[TURB_SOL]->OutputVariables[iPoint*nVar_Extra+jVar];
            }
            jPoint++;
            
          }
        }
        
        /*--- Gather the data on the master node. ---*/
        
#ifdef HAVE_MPI
        SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif
        
        /*--- The master node unpacks and sorts this variable by global index ---*/
        
        if (rank == MASTER_NODE) {
          jPoint = 0; iVar = iVar_Extra;
          for (iProcessor = 0; iProcessor < size; iProcessor++) {
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
  
  /*--- Release memory needed for surface coefficients ---*/
  
  delete [] Local_Halo;
  
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    delete [] Aux_Frict; delete [] Aux_Heat; delete [] Aux_yPlus;
  }
  if (( Kind_Solver == ADJ_EULER              ) ||
      ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
      ( Kind_Solver == ADJ_RANS               ) ||
      ( Kind_Solver == DISC_ADJ_EULER         ) ||
      ( Kind_Solver == DISC_ADJ_NAVIER_STOKES ) ||
      ( Kind_Solver == DISC_ADJ_RANS          )) {
    delete [] Aux_Sens;
  }
  
}

void COutput::MergeBaselineSolution(CConfig *config, CGeometry *geometry, CSolver *solver, unsigned short val_iZone) {
  
  /*--- Local variables needed on all processors ---*/
  unsigned short iVar;
  unsigned long iPoint = 0, jPoint = 0;
  
  nVar_Total = config->fields.size() - 1;
  
  /*--- Merge the solution either in serial or parallel. ---*/
  
#ifndef HAVE_MPI
  
  /*--- In serial, the single process has access to all solution data,
   so it is simple to retrieve and store inside Solution_Data. ---*/
  
  unsigned short iMarker;
  unsigned long iVertex, nTotalPoints = 0;
  int SendRecv;
  
  /*--- First, create a structure to locate any periodic halo nodes ---*/
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
            (SendRecv < 0)) {
          Local_Halo[iPoint] = false;
        }
      }
      
    }
  }
  
  /*--- Total number of points in the mesh (this might include periodic points). ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (!Local_Halo[iPoint]) nTotalPoints++;
  
  nGlobal_Poin = nTotalPoints;
  Data = new su2double*[nVar_Total];
  for (iVar = 0; iVar < nVar_Total; iVar++) {
    Data[iVar] = new su2double[nGlobal_Poin];
  }
  
  /*--- Loop over all points in the mesh, but only write data
   for nodes in the domain (ignore periodic halo nodes). ---*/
  
  jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    if (!Local_Halo[iPoint]) {
      
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
  
  int rank, nProcessor, iProcessor;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  /*--- Local variables needed for merging with MPI ---*/
  
  unsigned long iVertex, iMarker;
  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;
  
  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/
  
  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        
        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }
    
    /*--- Sum total number of nodes that belong to the domain ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;
    
  }
  Buffer_Send_nPoint[0] = nLocalPoint;
  
  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];
  
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  nBuffer_Scalar = MaxLocalPoint;
  
  /*--- Send and Recv buffers. ---*/
  
  su2double *Buffer_Send_Var = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Var = NULL;
  
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Prepare the receive buffers in the master node only. ---*/
  if (rank == MASTER_NODE) {
    
    Buffer_Recv_Var = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];
    
    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
    }
    Data = new su2double*[nVar_Total];
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      Data[iVar] = new su2double[nGlobal_Poin];
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
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
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

    SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    if (iVar == 0) {
      SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
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
  
  delete [] Local_Halo;
  
}

void COutput::SetRestart(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  /*--- Local variables ---*/
  
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iExtIter = config->GetExtIter();
  bool grid_movement = config->GetGrid_Movement();
  ofstream restart_file;
  string filename;
  
  /*--- Retrieve filename from config ---*/
  
  if ((config->GetAdjoint()) || (config->GetDiscrete_Adjoint())) {
    filename = config->GetRestart_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else {
    filename = config->GetRestart_FlowFileName();
    filename = config->GetRestart_FlowFileName(filename, val_iZone);
  }
  
  /*--- Unsteady problems require an iteration number to be appended. ---*/
  
  if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(val_iZone));
  } else if (config->GetWrt_Unsteady()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
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
	if ( Kind_Solver == LINEAR_ELASTICITY )
    restart_file << "\t\"Displacement_" << iVar+1<<"\"";
	else
    restart_file << "\t\"Conservative_" << iVar+1<<"\"";
  }

  if (!config->GetLow_MemoryOutput()) {
    
    if (config->GetWrt_Limiters()) {
      for (iVar = 0; iVar < nVar_Consv; iVar++) {
        restart_file << "\t\"Limiter_" << iVar+1<<"\"";
      }
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
      restart_file << "\t\"Pressure\"\t\"Temperature\"\t\"Pressure_Coefficient\"\t\"Mach\"";
    }
    
    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      restart_file << "\t\"Laminar_Viscosity\"\t\"Skin_Friction_Coefficient\"\t\"Heat_Flux\"\t\"Y_Plus\"";
    }
    
    if (Kind_Solver == RANS) {
      restart_file << "\t\"Eddy_Viscosity\"";
    }
    
    if (config->GetWrt_SharpEdges()) {
      if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
        restart_file << "\t\"Sharp_Edge_Dist\"";
      }
    }
    
    if (Kind_Solver == POISSON_EQUATION) {
      for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
        restart_file << "\t\"poissonField_" << iDim+1 << "\"";
    }
    
    if ((Kind_Solver == ADJ_EULER              ) ||
        (Kind_Solver == ADJ_NAVIER_STOKES      ) ||
        (Kind_Solver == ADJ_RANS               )   ) {
      restart_file << "\t\"Surface_Sensitivity\"\t\"Solution_Sensor\"";
    }
    if (( Kind_Solver == DISC_ADJ_EULER              ) ||
        ( Kind_Solver == DISC_ADJ_NAVIER_STOKES      ) ||
        ( Kind_Solver == DISC_ADJ_RANS               )) {
      restart_file << "\t\"Surface_Sensitivity\"\t\"Sensitivity_x\"\t\"Sensitivity_y\"";
      if (geometry->GetnDim() == 3){
        restart_file << "\t\"Sensitivity_z\"";
      }
    }
    
    if ((Kind_Solver == LINEAR_ELASTICITY) && (geometry->GetnDim() == 2)) {
      restart_file << "\t\"Sxx\"\t\"Syy\"\t\"Sxy\"\t\"Von_Mises_Stress\"\t\"Flow_Pressure\"";
    }

    if ((Kind_Solver == LINEAR_ELASTICITY) && (geometry->GetnDim() == 3)) {
      restart_file << "\t\"Sxx\"\t\"Syy\"\t\"Sxy\"\t\"Szz\"\t\"Sxz\"\t\"Syz\"\t\"Von_Mises_Stress\"\t\"Flow_Pressure\"";
    }
    
    if (config->GetExtraOutput()) {
      string *headings = NULL;
      //if (Kind_Solver == RANS) {
      headings = solver[TURB_SOL]->OutputHeadingNames;
      //}
      
      for (iVar = 0; iVar < nVar_Extra; iVar++) {
        if (headings == NULL) {
          restart_file << "\t\"ExtraOutput_" << iVar+1<<"\"";
        } else{
          restart_file << "\t\""<< headings[iVar] <<"\"";
        }
      }
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
  
  unsigned short iDim, nDim = geometry->GetnDim();

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
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
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
      if (nGlobal_Pris > 0) delete [] Conn_Pris;
      if (nGlobal_Pyra > 0) delete [] Conn_Pyra;
    }
    
  }
}

void COutput::DeallocateSolution(CConfig *config, CGeometry *geometry) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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

void COutput::SetConvHistory_Header(ofstream *ConvHist_file, CConfig *config) {
  char cstr[200], buffer[50], turb_resid[1000];
  unsigned short iMarker, iMarker_Monitoring;
  string Monitoring_Tag, monitoring_coeff, aeroelastic_coeff;
  
  bool rotating_frame = config->GetRotating_Frame();
  bool aeroelastic = config->GetAeroelastic_Simulation();
  bool equiv_area = config->GetEquivArea();
  bool turbulent = ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == ADJ_RANS) ||
                    (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool frozen_turb = config->GetFrozen_Visc();
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool inv_design = (config->GetInvDesign_Cp() || config->GetInvDesign_HeatFlux());
  bool output_1d = config->GetWrt_1D_Output();
  bool output_per_surface = false;
  bool output_massflow = (config->GetKind_ObjFunc() == MASS_FLOW_RATE);
  if (config->GetnMarker_Monitoring() > 1) output_per_surface = true;
  
  unsigned short direct_diff = config->GetDirectDiff();

  bool isothermal = false;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ))
      isothermal = true;
  
  /*--- Write file name with extension ---*/
  
  string filename = config->GetConv_FileName();
  strcpy (cstr, filename.data());
  
  if (config->GetWrt_Unsteady() && config->GetRestart()) {
    long iExtIter = config->GetUnst_RestartIter();
    if (SU2_TYPE::Int(iExtIter) < 10) SPRINTF (buffer, "_0000%d", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 10) && (SU2_TYPE::Int(iExtIter) < 100)) SPRINTF (buffer, "_000%d", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 100) && (SU2_TYPE::Int(iExtIter) < 1000)) SPRINTF (buffer, "_00%d", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d", SU2_TYPE::Int(iExtIter));
    if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d", SU2_TYPE::Int(iExtIter));
    strcat(cstr, buffer);
  }
  
  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == FIELDVIEW)) SPRINTF (buffer, ".dat");
  else if ((config->GetOutput_FileFormat() == TECPLOT_BINARY) ||
           (config->GetOutput_FileFormat() == FIELDVIEW_BINARY))  SPRINTF (buffer, ".plt");
  else if (config->GetOutput_FileFormat() == PARAVIEW)  SPRINTF (buffer, ".csv");
  strcat(cstr, buffer);
  
  ConvHist_file->open(cstr, ios::out);
  ConvHist_file->precision(15);
  
  /*--- Begin of the header ---*/
  
  char begin[]= "\"Iteration\"";
  
  /*--- Header for the coefficients ---*/
  
  char flow_coeff[]= ",\"CLift\",\"CDrag\",\"CSideForce\",\"CMx\",\"CMy\",\"CMz\",\"CFx\",\"CFy\",\"CFz\",\"CL/CD\"";
  char heat_coeff[]= ",\"HeatFlux_Total\",\"HeatFlux_Maximum\"";
  char equivalent_area_coeff[]= ",\"CEquivArea\",\"CNearFieldOF\"";
  char rotating_frame_coeff[]= ",\"CMerit\",\"CT\",\"CQ\"";
  char free_surface_coeff[]= ",\"CFreeSurface\"";
  char wave_coeff[]= ",\"CWave\"";
  char fea_coeff[]= ",\"CFEA\"";
  char adj_coeff[]= ",\"Sens_Geo\",\"Sens_Mach\",\"Sens_AoA\",\"Sens_Press\",\"Sens_Temp\",\"Sens_AoS\",\"Sens_BPress\"";
  char oneD_outputs[]= ",\"Avg_TotalPress\",\"Avg_Mach\",\"Avg_Temperature\",\"MassFlowRate\",\"FluxAvg_Pressure\",\"FluxAvg_Density\",\"FluxAvg_Velocity\",\"FluxAvg_Enthalpy\"";
  char Cp_inverse_design[]= ",\"Cp_Diff\"";
  char Heat_inverse_design[]= ",\"HeatFlux_Diff\"";
  char mass_flow_rate[] = ",\"MassFlowRate\"";
  char d_flow_coeff[] = ",\"D(CLift)\",\"D(CDrag)\",\"D(CSideForce)\",\"D(CMx)\",\"D(CMy)\",\"D(CMz)\",\"D(CFx)\",\"D(CFy)\",\"D(CFz)\",\"D(CL/CD)\"";
  
  /* Find the markers being monitored and create a header for them */
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
    monitoring_coeff += ",\"CLift_"  + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CDrag_"  + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CSideForce_" + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CL/CD_" + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CFx_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CFy_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CFz_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CMx_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CMy_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CMz_"    + Monitoring_Tag + "\"";
    aeroelastic_coeff += ",\"plunge_" + Monitoring_Tag + "\"";
    aeroelastic_coeff += ",\"pitch_"  + Monitoring_Tag + "\"";
  }
  
  /*--- Header for the residuals ---*/

  char flow_resid[]= ",\"Res_Flow[0]\",\"Res_Flow[1]\",\"Res_Flow[2]\",\"Res_Flow[3]\",\"Res_Flow[4]\"";
  char adj_flow_resid[]= ",\"Res_AdjFlow[0]\",\"Res_AdjFlow[1]\",\"Res_AdjFlow[2]\",\"Res_AdjFlow[3]\",\"Res_AdjFlow[4]\"";
  switch (config->GetKind_Turb_Model()) {
    case SA:	   SPRINTF (turb_resid, ",\"Res_Turb[0]\""); break;
    case SA_NEG: SPRINTF (turb_resid, ",\"Res_Turb[0]\""); break;
    case SST:   	SPRINTF (turb_resid, ",\"Res_Turb[0]\",\"Res_Turb[1]\""); break;
  }
  char adj_turb_resid[]= ",\"Res_AdjTurb[0]\"";
  char levelset_resid[]= ",\"Res_LevelSet\"";
  char adj_levelset_resid[]= ",\"Res_AdjLevelSet\"";
  char wave_resid[]= ",\"Res_Wave[0]\",\"Res_Wave[1]\"";
  char fea_resid[]= ",\"Res_FEA\"";
  char heat_resid[]= ",\"Res_Heat\"";
  
  /*--- End of the header ---*/
  
  char end[]= ",\"Linear_Solver_Iterations\",\"CFL_Number\",\"Time(min)\"\n";
  
  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == TECPLOT_BINARY) ||
      (config->GetOutput_FileFormat() == FIELDVIEW) ||
      (config->GetOutput_FileFormat() == FIELDVIEW_BINARY)) {
    ConvHist_file[0] << "TITLE = \"SU2 Simulation\"" << endl;
    ConvHist_file[0] << "VARIABLES = ";
  }
  
  /*--- Write the header, case depending ---*/
  switch (config->GetKind_Solver()) {
      
    case EULER : case NAVIER_STOKES: case RANS :
      ConvHist_file[0] << begin << flow_coeff;
      if (isothermal) ConvHist_file[0] << heat_coeff;
      if (equiv_area) ConvHist_file[0] << equivalent_area_coeff;
      if (inv_design) {
        ConvHist_file[0] << Cp_inverse_design;
        if (isothermal) ConvHist_file[0] << Heat_inverse_design;
      }
      if (rotating_frame) ConvHist_file[0] << rotating_frame_coeff;
      ConvHist_file[0] << flow_resid;
      if (turbulent) ConvHist_file[0] << turb_resid;
      if (aeroelastic) ConvHist_file[0] << aeroelastic_coeff;
      if (output_per_surface) ConvHist_file[0] << monitoring_coeff;
      if (output_1d) ConvHist_file[0] << oneD_outputs;
      if (output_massflow and !output_1d)  ConvHist_file[0]<< mass_flow_rate;
      if (direct_diff != NO_DERIVATIVE) ConvHist_file[0] << d_flow_coeff;
      ConvHist_file[0] << end;
      if (freesurface) {
        ConvHist_file[0] << begin << flow_coeff << free_surface_coeff;
        ConvHist_file[0] << flow_resid << levelset_resid << end;
      }

      break;
      
    case ADJ_EULER      : case ADJ_NAVIER_STOKES      : case ADJ_RANS:
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
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
  
  if (config->GetOutput_FileFormat() == TECPLOT ||
      config->GetOutput_FileFormat() == TECPLOT_BINARY ||
      config->GetOutput_FileFormat() == FIELDVIEW ||
      config->GetOutput_FileFormat() == FIELDVIEW_BINARY) {
    ConvHist_file[0] << "ZONE T= \"Convergence history\"" << endl;
  }
  
}


void COutput::SetConvHistory_Body(ofstream *ConvHist_file,
                                     CGeometry ***geometry,
                                     CSolver ****solver_container,
                                     CConfig **config,
                                     CIntegration ***integration,
                                     bool DualTime_Iteration,
                                     su2double timeused,
                                     unsigned short val_iZone) {
  
  bool output_1d  = config[val_iZone]->GetWrt_1D_Output();
  bool output_massflow = (config[val_iZone]->GetKind_ObjFunc() == MASS_FLOW_RATE);
  unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
  
  int rank;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  rank = MASTER_NODE;
#endif
  
  /*--- If 1-D outputs requested, calculated them. Requires info from all nodes,
   Get area-averaged and flux-averaged values at the specified surface ---*/
  
  if (output_1d) {
    switch (config[val_iZone]->GetKind_Solver()) {
      case EULER:                   case NAVIER_STOKES:                   case RANS:
      case ADJ_EULER:               case ADJ_NAVIER_STOKES:               case ADJ_RANS:
        OneDimensionalOutput(solver_container[val_iZone][FinestMesh][FLOW_SOL], geometry[val_iZone][FinestMesh], config[val_iZone]);
        break;
    }
  }
  if (output_massflow and !output_1d) {
    switch (config[val_iZone]->GetKind_Solver()) {
      case EULER:                   case NAVIER_STOKES:                   case RANS:
      case ADJ_EULER:               case ADJ_NAVIER_STOKES:               case ADJ_RANS:
        SetMassFlowRate(solver_container[val_iZone][FinestMesh][FLOW_SOL], geometry[val_iZone][FinestMesh], config[val_iZone]);
        break;
    }
  }

  /*--- Output using only the master node ---*/
  if (rank == MASTER_NODE) {
    
    unsigned long iIntIter = config[val_iZone]->GetIntIter();
    unsigned long iExtIter = config[val_iZone]->GetExtIter();
    
    /*--- WARNING: These buffers have hard-coded lengths. Note that you
     may have to adjust them to be larger if adding more entries. ---*/
    char begin[1000], direct_coeff[1000], surface_coeff[1000], aeroelastic_coeff[1000], monitoring_coeff[10000],
    adjoint_coeff[1000], flow_resid[1000], adj_flow_resid[1000], turb_resid[1000], trans_resid[1000],
    adj_turb_resid[1000], levelset_resid[1000], adj_levelset_resid[1000], wave_coeff[1000],
    heat_coeff[1000], fea_coeff[1000], wave_resid[1000], heat_resid[1000], fea_resid[1000], end[1000],
    oneD_outputs[1000], massflow_outputs[1000], d_direct_coeff[1000];

    su2double dummy = 0.0, *Coord;
    unsigned short iVar, iMarker, iMarker_Monitoring;
    
    unsigned long LinSolvIter = 0, iPointMaxResid;
    su2double timeiter = timeused/su2double(iExtIter+1);
    
    unsigned short nDim = geometry[val_iZone][FinestMesh]->GetnDim();
    
    bool compressible = (config[val_iZone]->GetKind_Regime() == COMPRESSIBLE);
    bool incompressible = (config[val_iZone]->GetKind_Regime() == INCOMPRESSIBLE);
    bool freesurface = (config[val_iZone]->GetKind_Regime() == FREESURFACE);
    
    bool rotating_frame = config[val_iZone]->GetRotating_Frame();
    bool aeroelastic = config[val_iZone]->GetAeroelastic_Simulation();
    bool equiv_area = config[val_iZone]->GetEquivArea();
    bool inv_design = (config[val_iZone]->GetInvDesign_Cp() || config[val_iZone]->GetInvDesign_HeatFlux());
    bool transition = (config[val_iZone]->GetKind_Trans_Model() == LM);
    bool isothermal = false;
    for (iMarker = 0; iMarker < config[val_iZone]->GetnMarker_All(); iMarker++)
      if ((config[val_iZone]->GetMarker_All_KindBC(iMarker) == ISOTHERMAL))
        isothermal = true;
    bool turbulent = ((config[val_iZone]->GetKind_Solver() == RANS) || (config[val_iZone]->GetKind_Solver() == ADJ_RANS) ||
                      (config[val_iZone]->GetKind_Solver() == DISC_ADJ_RANS));
    bool adjoint = config[val_iZone]->GetAdjoint() || config[val_iZone]->GetDiscrete_Adjoint();
    bool disc_adj = config[val_iZone]->GetDiscrete_Adjoint();
    bool wave = (config[val_iZone]->GetKind_Solver() == WAVE_EQUATION);
    bool heat = (config[val_iZone]->GetKind_Solver() == HEAT_EQUATION);
    bool fea = (config[val_iZone]->GetKind_Solver() == LINEAR_ELASTICITY);
    bool flow = (config[val_iZone]->GetKind_Solver() == EULER) || (config[val_iZone]->GetKind_Solver() == NAVIER_STOKES) ||
    (config[val_iZone]->GetKind_Solver() == RANS) || (config[val_iZone]->GetKind_Solver() == ADJ_EULER) ||
    (config[val_iZone]->GetKind_Solver() == ADJ_NAVIER_STOKES) || (config[val_iZone]->GetKind_Solver() == ADJ_RANS);
    
    bool turbo = config[val_iZone]->GetBoolTurboPerf();
    string inMarker_Tag, outMarker_Tag;

    bool output_per_surface = false;
    if (config[val_iZone]->GetnMarker_Monitoring() > 1) output_per_surface = true;

    unsigned short direct_diff = config[val_iZone]->GetDirectDiff();


    /*--- Initialize variables to store information from all domains (direct solution) ---*/
    su2double Total_CLift = 0.0, Total_CDrag = 0.0, Total_CSideForce = 0.0, Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0, Total_CEff = 0.0,
    Total_CEquivArea = 0.0, Total_CNearFieldOF = 0.0, Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0, Total_CMerit = 0.0,
    Total_CT = 0.0, Total_CQ = 0.0, Total_CFreeSurface = 0.0, Total_CWave = 0.0, Total_CHeat = 0.0, Total_CpDiff = 0.0, Total_HeatFluxDiff = 0.0,
    Total_CFEA = 0.0, Total_Heat = 0.0, Total_MaxHeat = 0.0, Total_Mdot = 0.0;
    su2double OneD_AvgStagPress = 0.0, OneD_AvgMach = 0.0, OneD_AvgTemp = 0.0, OneD_MassFlowRate = 0.0,
    OneD_FluxAvgPress = 0.0, OneD_FluxAvgDensity = 0.0, OneD_FluxAvgVelocity = 0.0, OneD_FluxAvgEntalpy = 0.0;
    
    /*--- Initialize variables to store information from all zone for turboperformance (direct solution) ---*/
    su2double *TotalStaticEfficiency = NULL,
    *TotalTotalEfficiency = NULL,
	*KineticEnergyLoss 	  = NULL,
	*TotalPressureLoss 	  = NULL,
	*MassFlowIn 		  = NULL,
	*MassFlowOut          = NULL,
	*FlowAngleIn          = NULL,
	*FlowAngleOut         = NULL,
	*EulerianWork         = NULL,
	*TotalEnthalpyIn      = NULL,
	*PressureRatio        = NULL,
	*PressureOut          = NULL,
	*EnthalpyOut          = NULL,
	*MachIn               = NULL,
	*MachOut              = NULL,
	*NormalMachIn         = NULL,
	*NormalMachOut        = NULL,
	*VelocityOutIs        = NULL;



    /*--- Initialize variables to store information from all domains (adjoint solution) ---*/
    su2double Total_Sens_Geo = 0.0, Total_Sens_Mach = 0.0, Total_Sens_AoA = 0.0;
    su2double Total_Sens_Press = 0.0, Total_Sens_Temp = 0.0, Total_Sens_BPress=0.0;
    
    /*--- Initialize variables to store information from all domains (direct differentiation) ---*/
    su2double D_Total_CLift = 0.0, D_Total_CDrag = 0.0, D_Total_CSideForce = 0.0, D_Total_CMx = 0.0, D_Total_CMy = 0.0, D_Total_CMz = 0.0, D_Total_CEff = 0.0, D_Total_CFx = 0.0, D_Total_CFy = 0.0, D_Total_CFz = 0.0;

    /*--- Residual arrays ---*/
    su2double *residual_flow         = NULL,
    *residual_turbulent    = NULL,
    *residual_transition   = NULL,
    *residual_levelset     = NULL;
    su2double *residual_adjflow      = NULL,
    *residual_adjturbulent = NULL,
    *residual_adjlevelset  = NULL;
    su2double *residual_wave         = NULL;
    su2double *residual_fea          = NULL;
    su2double *residual_heat         = NULL;
    
    /*--- Coefficients Monitored arrays ---*/
    su2double *aeroelastic_plunge = NULL,
    *aeroelastic_pitch  = NULL,
    *Surface_CLift      = NULL,
    *Surface_CDrag      = NULL,
    *Surface_CSideForce = NULL,
    *Surface_CEff       = NULL,
    *Surface_CFx        = NULL,
    *Surface_CFy        = NULL,
    *Surface_CFz        = NULL,
    *Surface_CMx        = NULL,
    *Surface_CMy        = NULL,
    *Surface_CMz        = NULL;
    
    /*--- Initialize number of variables ---*/
    unsigned short nVar_Flow = 0, nVar_LevelSet = 0, nVar_Turb = 0,
    nVar_Trans = 0, nVar_Wave = 0, nVar_Heat = 0, nVar_FEA = 0,
    nVar_AdjFlow = 0, nVar_AdjLevelSet = 0, nVar_AdjTurb = 0;
    
    /*--- Direct problem variables ---*/
    if (compressible) nVar_Flow = nDim+2; else nVar_Flow = nDim+1;
    if (turbulent) {
      switch (config[val_iZone]->GetKind_Turb_Model()) {
        case SA:	   nVar_Turb = 1; break;
        case SA_NEG: nVar_Turb = 1; break;
        case SST:    nVar_Turb = 2; break;
      }
    }
    if (transition) nVar_Trans = 2;
    if (wave) nVar_Wave = 2;
    if (fea) nVar_FEA = nDim;
    if (heat) nVar_Heat = 1;
    if (freesurface) nVar_LevelSet = 1;
    
    /*--- Adjoint problem variables ---*/
    if (compressible) nVar_AdjFlow = nDim+2; else nVar_AdjFlow = nDim+1;
    if (turbulent) {
      switch (config[val_iZone]->GetKind_Turb_Model()) {
        case SA:	   nVar_AdjTurb = 1; break;
        case SA_NEG: nVar_AdjTurb = 1; break;
        case SST:    nVar_AdjTurb = 2; break;
      }
    }
    if (freesurface) nVar_AdjLevelSet = 1;
    
    /*--- Allocate memory for the residual ---*/
    residual_flow       = new su2double[nVar_Flow];
    residual_turbulent  = new su2double[nVar_Turb];
    residual_transition = new su2double[nVar_Trans];
    residual_levelset   = new su2double[nVar_LevelSet];
    residual_wave       = new su2double[nVar_Wave];
    residual_fea        = new su2double[nVar_FEA];
    residual_heat       = new su2double[nVar_Heat];
    
    residual_adjflow      = new su2double[nVar_AdjFlow];
    residual_adjturbulent = new su2double[nVar_AdjTurb];
    residual_adjlevelset  = new su2double[nVar_AdjLevelSet];
    
    /*--- Allocate memory for the coefficients being monitored ---*/
    aeroelastic_plunge = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    aeroelastic_pitch  = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CLift      = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CDrag      = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CSideForce = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CEff       = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFx        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFy        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFz        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMx        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMy        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMz        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    
    /*--- Allocate memory for the turboperformace ---*/
    TotalStaticEfficiency = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    TotalTotalEfficiency  = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    KineticEnergyLoss 	  = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    TotalPressureLoss 	  = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    MassFlowIn 		      = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    MassFlowOut           = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    FlowAngleIn           = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    FlowAngleOut          = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    EulerianWork          = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    TotalEnthalpyIn       = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    PressureRatio         = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    PressureOut           = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    EnthalpyOut           = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    MachIn                = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    MachOut               = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    NormalMachIn          = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    NormalMachOut         = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];
    VelocityOutIs         = new su2double[config[ZONE_0]->Get_nMarkerTurboPerf()];





    /*--- Write information from nodes ---*/
    switch (config[val_iZone]->GetKind_Solver()) {
        
      case EULER:                   case NAVIER_STOKES:                   case RANS:
      case ADJ_EULER:               case ADJ_NAVIER_STOKES:               case ADJ_RANS:
      case DISC_ADJ_EULER:          case DISC_ADJ_NAVIER_STOKES:          case DISC_ADJ_RANS:
        
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

        if (direct_diff != NO_DERIVATIVE){
          D_Total_CLift       = SU2_TYPE::GetDerivative(Total_CLift);
          D_Total_CDrag       = SU2_TYPE::GetDerivative(Total_CDrag);
          D_Total_CSideForce  = SU2_TYPE::GetDerivative(Total_CSideForce);
          D_Total_CEff        = SU2_TYPE::GetDerivative(Total_CEff);
          D_Total_CMx         = SU2_TYPE::GetDerivative(Total_CMx);
          D_Total_CMy         = SU2_TYPE::GetDerivative(Total_CMy);
          D_Total_CMz         = SU2_TYPE::GetDerivative(Total_CMz);
          D_Total_CFx         = SU2_TYPE::GetDerivative(Total_CFx);
          D_Total_CFy         = SU2_TYPE::GetDerivative(Total_CFy);
          D_Total_CFz         = SU2_TYPE::GetDerivative(Total_CFz);
        }

        if (freesurface) {
          Total_CFreeSurface = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFreeSurface();
        }
        
        if (isothermal) {
          Total_Heat     = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_HeatFlux();
          Total_MaxHeat  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_MaxHeatFlux();
        }
        
        if (equiv_area) {
          Total_CEquivArea    = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CEquivArea();
          Total_CNearFieldOF  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CNearFieldOF();
          
          /*--- Note that there is a redefinition of the nearfield based functionals ---*/
          Total_CEquivArea    = config[val_iZone]->GetWeightCd()*Total_CDrag + (1.0-config[val_iZone]->GetWeightCd())*Total_CEquivArea;
          Total_CNearFieldOF  = config[val_iZone]->GetWeightCd()*Total_CDrag + (1.0-config[val_iZone]->GetWeightCd())*Total_CNearFieldOF;
        }
        
        if (inv_design) {
          Total_CpDiff  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CpDiff();
          if (isothermal) {
            Total_HeatFluxDiff = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_HeatFluxDiff();
          }
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
          }
        }
        
        if (output_per_surface) {
          /*--- Look over the markers being monitored and get the desired values ---*/
          for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
            Surface_CLift[iMarker_Monitoring]      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CLift(iMarker_Monitoring);
            Surface_CDrag[iMarker_Monitoring]      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CDrag(iMarker_Monitoring);
            Surface_CSideForce[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CSideForce(iMarker_Monitoring);
            Surface_CEff[iMarker_Monitoring]       = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CEff(iMarker_Monitoring);
            Surface_CFx[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFx(iMarker_Monitoring);
            Surface_CFy[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFy(iMarker_Monitoring);
            Surface_CFz[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFz(iMarker_Monitoring);
            Surface_CMx[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMx(iMarker_Monitoring);
            Surface_CMy[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMy(iMarker_Monitoring);
            Surface_CMz[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMz(iMarker_Monitoring);
          }
        }
        
        if (turbo) {
        	/*--- Loop over the nMarker of turboperformance and get the desired values ---*/
        	for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->Get_nMarkerTurboPerf(); iMarker_Monitoring++) {
        		TotalStaticEfficiency[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotalStaticEfficiency(iMarker_Monitoring);
						TotalTotalEfficiency[iMarker_Monitoring]  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotalTotalEfficiency(iMarker_Monitoring);
						KineticEnergyLoss[iMarker_Monitoring] 	  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetKineticEnergyLoss(iMarker_Monitoring);
						TotalPressureLoss[iMarker_Monitoring] 	  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotalPressureLoss(iMarker_Monitoring);
						MassFlowIn[iMarker_Monitoring] 		      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMassFlowIn(iMarker_Monitoring);
						MassFlowOut[iMarker_Monitoring]           = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMassFlowOut(iMarker_Monitoring);
						FlowAngleIn[iMarker_Monitoring]           = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetFlowAngleIn(iMarker_Monitoring);
						FlowAngleOut[iMarker_Monitoring]          = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetFlowAngleOut(iMarker_Monitoring);
						EulerianWork[iMarker_Monitoring]          = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetEulerianWork(iMarker_Monitoring);
						TotalEnthalpyIn[iMarker_Monitoring]       = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotalEnthalpyIn(iMarker_Monitoring);
						PressureRatio[iMarker_Monitoring]         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPressureRatio(iMarker_Monitoring);
						PressureOut[iMarker_Monitoring]           = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPressureOut(iMarker_Monitoring);
						EnthalpyOut[iMarker_Monitoring]           = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetEnthalpyOut(iMarker_Monitoring);
						MachIn[iMarker_Monitoring]                = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMachIn(iMarker_Monitoring);
						MachOut[iMarker_Monitoring]               = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMachOut(iMarker_Monitoring);
						NormalMachIn[iMarker_Monitoring]          = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetNormalMachIn(iMarker_Monitoring);
						NormalMachOut[iMarker_Monitoring]         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetNormalMachOut(iMarker_Monitoring);
						VelocityOutIs[iMarker_Monitoring]         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetVelocityOutIs(iMarker_Monitoring);
        	}
        }


//        if (fluid_structure) {
//          Total_CFEA  = solver_container[ZONE_0][FinestMesh][FEA_SOL]->GetTotal_CFEA();
//        }
        
        if (output_1d) {
          
          /*--- Get area-averaged and flux-averaged values at the specified surface ---*/
          
          OneD_AvgStagPress = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_TotalPress();
          OneD_AvgMach = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_Mach();
          OneD_AvgTemp = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_Temp();
          OneD_MassFlowRate = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_MassFlowRate();
          
          OneD_FluxAvgPress = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_FluxAvgPress();
          OneD_FluxAvgDensity = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_FluxAvgDensity();
          OneD_FluxAvgVelocity = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_FluxAvgVelocity();
          OneD_FluxAvgEntalpy = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_FluxAvgEntalpy();
          
        }
        /*--- Get Mass Flow at the Monitored Markers ---*/


        if (output_massflow) {
          Total_Mdot = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_MassFlowRate();
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
//        if (fluid_structure) {
//          for (iVar = 0; iVar < nVar_FEA; iVar++)
//            residual_fea[iVar] = solver_container[ZONE_0][FinestMesh][FEA_SOL]->GetRes_RMS(iVar);
//        }
        
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
          Total_Sens_BPress  = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_BPress();
          
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
    
    /*--- Header frequency ---*/
    
    bool Unsteady = ((config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                     (config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND));
    bool In_NoDualTime = (!DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    bool In_DualTime_0 = (DualTime_Iteration && (iIntIter % config[val_iZone]->GetWrt_Con_Freq_DualTime() == 0));
    bool In_DualTime_1 = (!DualTime_Iteration && Unsteady);
    bool In_DualTime_2 = (Unsteady && DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    bool In_DualTime_3 = (Unsteady && !DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    
    bool write_heads;
    if (Unsteady) write_heads = (iIntIter == 0);
    else write_heads = (((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*40)) == 0));
    bool write_turbo = (((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*200)) == 0));
    if ((In_NoDualTime || In_DualTime_0 || In_DualTime_1) && (In_NoDualTime || In_DualTime_2 || In_DualTime_3)) {
      
      /*--- Prepare the history file output, note that the dual
       time output don't write to the history file ---*/
      if (!DualTime_Iteration) {
        
        /*--- Write the begining of the history file ---*/
        SPRINTF (begin, "%12d", SU2_TYPE::Int(iExtIter));
        
        /*--- Write the end of the history file ---*/
        SPRINTF (end, ", %12.10f, %12.10f, %12.10f\n", su2double(LinSolvIter), config[val_iZone]->GetCFL(MESH_0), timeused/60.0);
        
        /*--- Write the solution and residual of the history file ---*/
        switch (config[val_iZone]->GetKind_Solver()) {
            
          case EULER : case NAVIER_STOKES: case RANS:
          case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS: case DISC_ADJ_EULER:
          case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
            
            /*--- Direct coefficients ---*/
            SPRINTF (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                     Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
                     Total_CFz, Total_CEff);
            if (direct_diff != NO_DERIVATIVE){
              SPRINTF (d_direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                       D_Total_CLift, D_Total_CDrag, D_Total_CSideForce, D_Total_CMx, D_Total_CMy, D_Total_CMz, D_Total_CFx, D_Total_CFy,
                       D_Total_CFz, D_Total_CEff);
            }
            if (isothermal)
              SPRINTF (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy,
                       Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_Heat, Total_MaxHeat);
            if (equiv_area)
              SPRINTF (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CEquivArea, Total_CNearFieldOF);
            if (inv_design) {
              SPRINTF (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CpDiff);
              Total_CpDiff  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CpDiff();
              if (isothermal) {
                SPRINTF (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_Heat, Total_MaxHeat, Total_CpDiff, Total_HeatFluxDiff);
              }
            }
            if (rotating_frame)
              SPRINTF (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
                       Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CMerit, Total_CT, Total_CQ);
            
            if (freesurface) {
              SPRINTF (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
                       Total_CFz, Total_CEff, Total_CFreeSurface);
            }
//            if (fluid_structure)
//              SPRINTF (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz,
//                       Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CFEA);
            
            if (aeroelastic) {
              for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
                //Append one by one the surface coeff to aeroelastic coeff. (Think better way do this, maybe use string)
                if (iMarker_Monitoring == 0) {
                  SPRINTF(aeroelastic_coeff, ", %12.10f", aeroelastic_plunge[iMarker_Monitoring]);
                }
                else {
                  SPRINTF(surface_coeff, ", %12.10f", aeroelastic_plunge[iMarker_Monitoring]);
                  strcat(aeroelastic_coeff, surface_coeff);
                }
                SPRINTF(surface_coeff, ", %12.10f", aeroelastic_pitch[iMarker_Monitoring]);
                strcat(aeroelastic_coeff, surface_coeff);
              }
            }
            
            if (output_per_surface) {
              for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
                //Append one by one the surface coeff to monitoring coeff. (Think better way do this, maybe use string)
                if (iMarker_Monitoring == 0) {
                  SPRINTF(monitoring_coeff, ", %12.10f", Surface_CLift[iMarker_Monitoring]);
                }
                else {
                  SPRINTF(surface_coeff, ", %12.10f", Surface_CLift[iMarker_Monitoring]);
                  strcat(monitoring_coeff, surface_coeff);
                }
                SPRINTF(surface_coeff, ", %12.10f", Surface_CDrag[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CSideForce[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CEff[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CFx[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CFy[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CFz[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CMx[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CMy[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CMz[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
              }
            }
            
            
            /*--- Flow residual ---*/
            if (nDim == 2) {
              if (compressible) SPRINTF (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), dummy);
              if (incompressible || freesurface) SPRINTF (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), dummy, dummy);
            }
            else {
              if (compressible) SPRINTF (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), log10 (residual_flow[4]) );
              if (incompressible || freesurface) SPRINTF (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), dummy);
            }
            
            /*--- Turbulent residual ---*/
            if (turbulent) {
              switch(nVar_Turb) {
                case 1: SPRINTF (turb_resid, ", %12.10f", log10 (residual_turbulent[0])); break;
                case 2: SPRINTF (turb_resid, ", %12.10f, %12.10f", log10(residual_turbulent[0]), log10(residual_turbulent[1])); break;
              }
            }
            /*---- Averaged stagnation pressure at an exit ----*/
            if (output_1d) {
              SPRINTF( oneD_outputs, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", OneD_AvgStagPress, OneD_AvgMach, OneD_AvgTemp, OneD_MassFlowRate, OneD_FluxAvgPress, OneD_FluxAvgDensity, OneD_FluxAvgVelocity, OneD_FluxAvgEntalpy);
            }
            if (output_massflow && !output_1d) {
              SPRINTF(massflow_outputs,", %12.10f", Total_Mdot);
            }

            
            /*--- Transition residual ---*/
            if (transition) {
              SPRINTF (trans_resid, ", %12.10f, %12.10f", log10(residual_transition[0]), log10(residual_transition[1]));
            }
            
            /*--- Free surface residual ---*/
            if (freesurface) {
              SPRINTF (levelset_resid, ", %12.10f", log10 (residual_levelset[0]));
            }
            
            /*--- Fluid structure residual ---*/
//            if (fluid_structure) {
//              if (nDim == 2) SPRINTF (levelset_resid, ", %12.10f, %12.10f, 0.0", log10 (residual_fea[0]), log10 (residual_fea[1]));
//              else SPRINTF (levelset_resid, ", %12.10f, %12.10f, %12.10f", log10 (residual_fea[0]), log10 (residual_fea[1]), log10 (residual_fea[2]));
//            }
            
            if (adjoint) {
              
              /*--- Adjoint coefficients ---*/
              SPRINTF (adjoint_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, 0.0, %12.10f", Total_Sens_Geo, Total_Sens_Mach, Total_Sens_AoA, Total_Sens_Press, Total_Sens_Temp, Total_Sens_BPress);
              
              /*--- Adjoint flow residuals ---*/
              if (nDim == 2) {
                if (compressible) SPRINTF (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0", log10 (residual_adjflow[0]), log10 (residual_adjflow[1]), log10 (residual_adjflow[2]), log10 (residual_adjflow[3]) );
                if (incompressible || freesurface) SPRINTF (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, 0.0, 0.0", log10 (residual_adjflow[0]), log10 (residual_adjflow[1]), log10 (residual_adjflow[2]) );
              }
              else {
                if (compressible) SPRINTF (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_adjflow[0]), log10 (residual_adjflow[1]), log10 (residual_adjflow[2]), log10 (residual_adjflow[3]), log10 (residual_adjflow[4]) );
                if (incompressible || freesurface) SPRINTF (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0", log10 (residual_adjflow[0]), log10 (residual_adjflow[1]), log10 (residual_adjflow[2]), log10 (residual_adjflow[3]) );
              }
              
              /*--- Adjoint turbulent residuals ---*/
              if (turbulent)
                if (!config[val_iZone]->GetFrozen_Visc())
                  SPRINTF (adj_turb_resid, ", %12.10f", log10 (residual_adjturbulent[0]));
              
              /*--- Adjoint free surface residuals ---*/
              if (freesurface) SPRINTF (adj_levelset_resid, ", %12.10f", log10 (residual_adjlevelset[0]));
            }
            
            break;
            
          case WAVE_EQUATION:
            
            SPRINTF (direct_coeff, ", %12.10f", Total_CWave);
            SPRINTF (wave_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_wave[0]), log10 (residual_wave[1]), dummy, dummy, dummy );
            
            break;
            
          case HEAT_EQUATION:
            
            SPRINTF (direct_coeff, ", %12.10f", Total_CHeat);
            SPRINTF (heat_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_heat[0]), dummy, dummy, dummy, dummy );
            
            break;
            
          case LINEAR_ELASTICITY:
            
            SPRINTF (direct_coeff, ", %12.10f", Total_CFEA);
            SPRINTF (fea_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_fea[0]), dummy, dummy, dummy, dummy );
            
            break;
            
        }
      }
      
      /*--- Write the screen header---*/
      if ((write_heads) && !(!DualTime_Iteration && Unsteady)) {
        
        if (!Unsteady) {
          switch (config[val_iZone]->GetKind_Solver()) {
            case EULER : case NAVIER_STOKES: case RANS:
            case ADJ_EULER : case ADJ_NAVIER_STOKES: case ADJ_RANS:
              
              cout << endl << "---------------------- Local Time Stepping Summary ----------------------" << endl;
              
              for (unsigned short iMesh = FinestMesh; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++)
                cout << "MG level: "<< iMesh << " -> Min. DT: " << solver_container[val_iZone][iMesh][FLOW_SOL]->GetMin_Delta_Time()<<
                ". Max. DT: " << solver_container[val_iZone][iMesh][FLOW_SOL]->GetMax_Delta_Time() <<
                ". CFL: " << config[val_iZone]->GetCFL(iMesh)  << "." << endl;
              
              cout << "-------------------------------------------------------------------------" << endl;

              if (direct_diff != NO_DERIVATIVE){
                cout << endl << "---------------------- Direct Differentiation Summary -------------------" << endl;
                cout << "Coefficients are differentiated with respect to ";
                switch (direct_diff) {
                  case D_MACH:
                    cout << "Mach number." << endl;
                    break;
                  case D_AOA:
                    cout << "AoA." << endl;
                    break;
                  case D_SIDESLIP:
                    cout << "AoS." << endl;
                    break;
                  case D_REYNOLDS:
                    cout << "Reynolds number." << endl;
                    break;
                  case D_TURB2LAM:
                    cout << "Turb/Lam ratio." << endl;
                    break;
                  case D_PRESSURE:
                    cout << "Freestream Pressure." << endl;
                    break;
                  case D_TEMPERATURE:
                    cout << "Freestream Temperature." << endl;
                    break;
                  case D_DENSITY:
                    cout << "Freestream Density." << endl;
                    break;
                  case D_VISCOSITY:
                    cout << "Freestream Viscosity." << endl;
                    break;
                  case D_DESIGN:
                    cout << "Design Variables." << endl;
                    break;
                  default:
                    break;
                  }

                cout << "    D_CLift(Total)" << "    D_CDrag(Total)" << "      D_CMz(Total)" <<"     D_CEff(Total)" << endl;
                cout.width(18); cout << D_Total_CLift;
                cout.width(18); cout << D_Total_CDrag;
                cout.width(18); cout << D_Total_CMz;
                cout.width(18); cout << D_Total_CEff;
                cout << endl << "-------------------------------------------------------------------------" << endl;
                cout << endl;
              }
              if (turbo && write_turbo){
              	cout << endl << "---------------------- Turbo Performance Summary -------------------" << endl;
              	for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->Get_nMarkerTurboPerf(); iMarker_Monitoring++){
              		inMarker_Tag = config[ZONE_0]->GetMarker_TurboPerf_BoundIn(iMarker_Monitoring);
              		outMarker_Tag = config[ZONE_0]->GetMarker_TurboPerf_BoundOut(iMarker_Monitoring);
              		switch (config[ZONE_0]->GetKind_TurboPerf(iMarker_Monitoring)) {
              			case BLADE:
											cout << "Blade performance between boundaries " << inMarker_Tag << " and "<< outMarker_Tag << " : "<<endl;
											cout << endl;
											cout << "   Total Pressure Loss(%)" << "   Kinetic Energy Loss(%)" << "            Eulerian Work" << endl;
											cout.width(25); cout << TotalPressureLoss[iMarker_Monitoring]*100.0;
											cout.width(25); cout << KineticEnergyLoss[iMarker_Monitoring]*100.0;
											cout.width(25); cout << EulerianWork[iMarker_Monitoring]*config[ZONE_0]->GetEnergy_Ref();
											cout << endl;
											cout << endl;
											cout << "     Total Inlet Enthalpy" << "          Outlet Enthalpy" << "            D_MassFlow(%)" <<  endl;
											cout.width(25); cout << TotalEnthalpyIn[iMarker_Monitoring]*config[ZONE_0]->GetEnergy_Ref();
											cout.width(25); cout << EnthalpyOut[iMarker_Monitoring]*config[ZONE_0]->GetEnergy_Ref();
											cout.width(25); cout << abs((MassFlowIn[iMarker_Monitoring] + MassFlowOut[iMarker_Monitoring])/MassFlowIn[iMarker_Monitoring])*100.0;
											cout << endl;
											cout << endl;
											cout << "   Isentropic Outlet Vel." << "         Inlet Flow Angle" << "        Outlet Flow Angle" <<endl;
											cout.width(25); cout << VelocityOutIs[iMarker_Monitoring]*config[ZONE_0]->GetVelocity_Ref();
											cout.width(25); cout << 180.0/PI_NUMBER*FlowAngleIn[iMarker_Monitoring];
											cout.width(25); cout << 180.0/PI_NUMBER*FlowAngleOut[iMarker_Monitoring];
											cout << endl;
											cout << endl;
											cout << "          Inlet Mass Flow"<< "               Inlet Mach" << "              Outlet Mach" << endl;
											cout.width(25); cout << MassFlowIn[iMarker_Monitoring]*config[ZONE_0]->GetVelocity_Ref()*config[ZONE_0]->GetDensity_Ref();
											cout.width(25); cout << MachIn[iMarker_Monitoring];
											cout.width(25); cout << MachOut[iMarker_Monitoring];
											cout << endl;
											cout << endl;
											cout << "        Inlet Normal Mach" << "       Outlet Normal Mach" << endl;
											cout.width(25); cout << NormalMachIn[iMarker_Monitoring];
											cout.width(25); cout << NormalMachOut[iMarker_Monitoring];
											cout << endl;
											cout << endl;
											cout << "           Pressure Ratio" << "         Outlet Pressure" << endl;
											cout.width(25); cout << PressureRatio[iMarker_Monitoring];
											cout.width(25); cout << PressureOut[iMarker_Monitoring]*config[ZONE_0]->GetPressure_Ref();
											cout << endl;
											cout << endl << "-------------------------------------------------------------------------" << endl;
											cout << endl;

											break;
              			case STAGE:
											cout << "Stage performance between boundaries " << inMarker_Tag << " and "<< outMarker_Tag << " : "<<endl;
											cout << endl;
											cout << "    Tot Tot Efficiency(%)" << "   Tot Stat Efficiency(%)" << endl;
											cout.width(25); cout << TotalTotalEfficiency[iMarker_Monitoring]*100.0;
											cout.width(25); cout << TotalStaticEfficiency[iMarker_Monitoring]*100.0;
											cout << endl;
											cout << endl;
											cout << "           Pressure Ratio" << "          Outlet Pressure" << endl;
											cout.width(25); cout << PressureRatio[iMarker_Monitoring];
											cout.width(25); cout << PressureOut[iMarker_Monitoring];
											cout << endl;
											cout << endl;
											cout << "     Total Inlet Enthalpy" << "    Total Outlet Enthalpy" << endl;
											cout.width(25); cout << TotalEnthalpyIn[iMarker_Monitoring];
											cout.width(25); cout << EnthalpyOut[iMarker_Monitoring];
											cout << endl;
											cout << endl << "-------------------------------------------------------------------------" << endl;
											cout << endl;

											break;
              			case TURBINE:
											cout << "Multi-stage performance between boundaries " << inMarker_Tag << " and "<< outMarker_Tag << " : "<<endl;
											cout << endl;
											cout << "    Tot Tot Efficiency(%)" << "   Tot Stat Efficiency(%)" << endl;
											cout.width(25); cout << TotalTotalEfficiency[iMarker_Monitoring]*100.0;
											cout.width(25); cout << TotalStaticEfficiency[iMarker_Monitoring]*100.0;
											cout << endl;
											cout << endl;
											cout << "           Pressure Ratio" << "          Outlet Pressure" << endl;
											cout.width(25); cout << PressureRatio[iMarker_Monitoring];
											cout.width(25); cout << PressureOut[iMarker_Monitoring];
											cout << endl;
											cout << endl;
											cout << "     Total Inlet Enthalpy" << "    Total Outlet Enthalpy" << endl;
											cout.width(25); cout << TotalEnthalpyIn[iMarker_Monitoring];
											cout.width(25); cout << EnthalpyOut[iMarker_Monitoring];
											cout << endl;
											cout << endl << "-------------------------------------------------------------------------" << endl;
											cout << endl;
											break;
              			default:
              				break;
              		}
              	}


              }
              break;

            case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
               cout << endl;
               cout << "------------------------ Discrete Adjoint Summary -----------------------" << endl;
               cout << "Total Geometry Sensitivity (updated every "  << config[val_iZone]->GetWrt_Sol_Freq() << " iterations): ";
               cout.precision(4);
               cout.setf(ios::scientific, ios::floatfield);
               cout << Total_Sens_Geo;
               cout << endl << "-------------------------------------------------------------------------" << endl;
              break;

          }
        }
        else {
          if (flow) {
            cout << endl << "Min DT: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<<
            ".Max DT: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() <<
            ".Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
          }
          else {
            cout << endl << "Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
          }
        }
        
        switch (config[val_iZone]->GetKind_Solver()) {
          case EULER :                  case NAVIER_STOKES:
            
            /*--- Visualize the maximum residual ---*/
            iPointMaxResid = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max(0);
            Coord = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max_Coord(0);
            
            cout << endl << "----------------------- Residual Evolution Summary ----------------------" << endl;

            cout << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_Max(0)) << "." << endl;
            
            if (config[val_iZone]->GetSystemMeasurements() == SI) {
              cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] << ", " << Coord[1];
              if (nDim == 3) cout << ", " << Coord[2];
              cout <<   ")." << endl;
            }
            else {
              cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0]*12.0 << ", " << Coord[1]*12.0;
              if (nDim == 3) cout << ", " << Coord[2]*12.0;
              cout <<   ")." << endl;
            }
            
            /*--- Print out the number of non-physical points and reconstructions ---*/
            
            if (config[val_iZone]->GetNonphysical_Points() > 0)
              cout << "There are " << config[val_iZone]->GetNonphysical_Points() << " non-physical points in the solution." << endl;
            if (config[val_iZone]->GetNonphysical_Reconstr() > 0)
              cout << "There are " << config[val_iZone]->GetNonphysical_Reconstr() << " non-physical states in the upwind reconstruction." << endl;
            
            cout << "-------------------------------------------------------------------------" << endl;

            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            
//            if (!fluid_structure) {
              if (incompressible) cout << "   Res[Press]" << "     Res[Velx]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
              else if (freesurface) cout << "   Res[Press]" << "     Res[Dist]" << "   CLift(Total)" << "     CLevelSet" << endl;
              else if (rotating_frame && nDim == 3 && !turbo) cout << "     Res[Rho]" << "     Res[RhoE]" << " CThrust(Total)" << " CTorque(Total)" << endl;
              else if (aeroelastic) cout << "     Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << "         plunge" << "          pitch" << endl;
              else if (equiv_area) cout << "     Res[Rho]" << "   CLift(Total)" << "   CDrag(Total)" << "    CPress(N-F)" << endl;
              else if (turbo)
            	  switch (config[ZONE_0]->GetKind_TurboPerf(0)) {
            	  	case BLADE:
            	  		cout << "     Res[Rho]" << "     Res[RhoE]"  << "  KineticLoss(%)" << "  D_MassFlow(%)" << endl;
            	  		break;
            	  	case STAGE: case TURBINE:
            	  		cout << "     Res[Rho]" << "     Res[RhoE]"  << " TSEfficiency(%)" << " Outlet Pressure" << endl;
            	  		break;
            	  	default:
            	  		break;
            	  }
              else cout << "     Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
//            }
//            else if (fluid_structure) cout << "     Res[Rho]" << "   Res[Displx]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
            
            break;
            
          case RANS :
            
            /*--- Visualize the maximum residual ---*/
            iPointMaxResid = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max(0);
            Coord = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max_Coord(0);
            
            cout << endl << "----------------------- Residual Evolution Summary ----------------------" << endl;

            cout << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_Max(0)) << "." << endl;
            if (config[val_iZone]->GetSystemMeasurements() == SI) {
              cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] << ", " << Coord[1];
              if (nDim == 3) cout << ", " << Coord[2];
              cout <<   ")." << endl;
            }
            else {
              cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0]*12.0 << ", " << Coord[1]*12.0;
              if (nDim == 3) cout << ", " << Coord[2]*12.0;
              cout <<   ")." << endl;
            }
            cout <<"Maximum Omega " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOmega_Max() << ", maximum Strain Rate " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetStrainMag_Max() << "." << endl;
            
            /*--- Print out the number of non-physical points and reconstructions ---*/
            if (config[val_iZone]->GetNonphysical_Points() > 0)
              cout << "There are " << config[val_iZone]->GetNonphysical_Points() << " non-physical points in the solution." << endl;
            if (config[val_iZone]->GetNonphysical_Reconstr() > 0)
              cout << "There are " << config[val_iZone]->GetNonphysical_Reconstr() << " non-physical states in the upwind reconstruction." << endl;
            
            cout << "-------------------------------------------------------------------------" << endl;

            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            if (incompressible || freesurface) cout << "   Res[Press]";
            else cout << "      Res[Rho]";//, cout << "     Res[RhoE]";
            
            switch (config[val_iZone]->GetKind_Turb_Model()) {
              case SA:	   cout << "       Res[nu]"; break;
              case SA_NEG: cout << "       Res[nu]"; break;
              case SST:	   cout << "     Res[kine]" << "     Res[omega]"; break;
            }
            
            if (transition) { cout << "      Res[Int]" << "       Res[Re]"; }
            else if (rotating_frame && nDim == 3 && !turbo ) cout << "   CThrust(Total)" << "   CTorque(Total)" << endl;
            else if (aeroelastic) cout << "   CLift(Total)" << "   CDrag(Total)" << "         plunge" << "          pitch" << endl;
            else if (equiv_area) cout << "   CLift(Total)" << "   CDrag(Total)" << "    CPress(N-F)" << endl;
            else if (turbo)
            	switch (config[ZONE_0]->GetKind_TurboPerf(0)) {
            		case BLADE:
            			cout << "  KineticLoss(%)" << "  D_MassFlow(%)" << endl;
									break;
            		case STAGE: case TURBINE:
            			cout << " TSEfficiency(%)" << " Outlet Pressure" << endl;
            			break;
            	  default:
									break;
            	}

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
          case DISC_ADJ_EULER:          case DISC_ADJ_NAVIER_STOKES:
            
            /*--- Visualize the maximum residual ---*/
            iPointMaxResid = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max(0);
            Coord = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max_Coord(0);
            cout << endl << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_Max(0)) << "." << endl;
            if (config[val_iZone]->GetSystemMeasurements() == SI) {
              cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] << ", " << Coord[1];
              if (nDim == 3) cout << ", " << Coord[2];
              cout <<   ")." << endl;
            }
            else {
              cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0]*12.0 << ", " << Coord[1]*12.0;
              if (nDim == 3) cout << ", " << Coord[2]*12.0;
              cout <<   ")." << endl;
            }
            
            /*--- Print out the number of non-physical points and reconstructions ---*/
            if (config[val_iZone]->GetNonphysical_Points() > 0)
              cout << "There are " << config[val_iZone]->GetNonphysical_Points() << " non-physical points in the solution." << endl;

            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            if (incompressible || freesurface) cout << "   Res[Psi_Press]" << "   Res[Psi_Velx]";
            else cout << "   Res[Psi_Rho]" << "     Res[Psi_E]";
            if (disc_adj){
              cout << "    Sens_Press" << "     Sens_Mach" << endl;
            } else {
              if (output_1d)
                cout << "      Sens_Geo" << "     Sens_BPress" << endl;
              else
                cout << "      Sens_Geo" << "     Sens_Mach" << endl;
            }
            if (freesurface) {
              cout << "   Res[Psi_Press]" << "   Res[Psi_Dist]" << "    Sens_Geo";
              if (output_1d)
                cout << "     Sens_BPress" << endl;
              else
                cout << "   Sens_Mach" << endl;
            }
            break;
            
          case ADJ_RANS : case DISC_ADJ_RANS:
            
            /*--- Visualize the maximum residual ---*/
            iPointMaxResid = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max(0);
            Coord = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max_Coord(0);
            cout << endl << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_Max(0)) << "." << endl;
            if (config[val_iZone]->GetSystemMeasurements() == SI) {
              cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] << ", " << Coord[1];
              if (nDim == 3) cout << ", " << Coord[2];
              cout <<   ")." << endl;
            }
            else {
              cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0]*12.0 << ", " << Coord[1]*12.0;
              if (nDim == 3) cout << ", " << Coord[2]*12.0;
              cout <<   ")." << endl;
            }
            
            /*--- Print out the number of non-physical points and reconstructions ---*/
            if (config[val_iZone]->GetNonphysical_Points() > 0)
              cout << "There are " << config[val_iZone]->GetNonphysical_Points() << " non-physical points in the solution." << endl;

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
            if (disc_adj){
              cout << "    Sens_Press" << "     Sens_Mach" << endl;
            } else {
              if (output_1d)
               cout << "      Sens_Geo" << "     Sens_BPress" << endl;
             else
               cout << "      Sens_Geo" << "     Sens_Mach" << endl;
            }
            if (freesurface) {
              cout << "   Res[Psi_Press]" << "   Res[Psi_Dist]" << "    Sens_Geo";
              if (output_1d)
                cout << "     Sens_BPress" << endl;
              else
                cout << "   Sens_Mach" << endl;
            }
            break;
            
        }
        
      }
      
      /*--- Write the solution on the screen and history file ---*/
      cout.precision(6);
      cout.setf(ios::fixed, ios::floatfield);
      
      if (!Unsteady) {
        cout.width(5); cout << iExtIter;
        cout.width(11); cout << timeiter;
        
      } else {
        cout.width(8); cout << iIntIter;
        cout.width(8); cout << iExtIter;
      }
      
      switch (config[val_iZone]->GetKind_Solver()) {
        case EULER : case NAVIER_STOKES:
          
          if (!DualTime_Iteration) {
            if (compressible) ConvHist_file[0] << begin << direct_coeff << flow_resid;
            if (incompressible) ConvHist_file[0] << begin << direct_coeff << flow_resid;
            if (freesurface) ConvHist_file[0] << begin << direct_coeff << flow_resid << levelset_resid << end;
//            if (fluid_structure) ConvHist_file[0] << fea_resid;
            if (aeroelastic) ConvHist_file[0] << aeroelastic_coeff;
            if (output_per_surface) ConvHist_file[0] << monitoring_coeff;
            if (output_1d) ConvHist_file[0] << oneD_outputs;
            if (output_massflow and !output_1d) ConvHist_file[0] << massflow_outputs;
            if (direct_diff != NO_DERIVATIVE) ConvHist_file[0] << d_direct_coeff;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed, ios::floatfield);
          cout.width(13); cout << log10(residual_flow[0]);
//          if (!fluid_structure && !equiv_area) {
          if (!equiv_area) {
            if (compressible) {
              if (nDim == 2 ) { cout.width(14); cout << log10(residual_flow[3]); }
              else { cout.width(14); cout << log10(residual_flow[4]); }
            }
            if (incompressible) { cout.width(14); cout << log10(residual_flow[1]); }
            if (freesurface) { cout.width(14); cout << log10(residual_levelset[0]); }
          }
//          else if (fluid_structure) { cout.width(14); cout << log10(residual_fea[0]); }
          
          if (rotating_frame && nDim == 3 && !turbo ) {
            cout.setf(ios::scientific, ios::floatfield);
            cout.width(15); cout << Total_CT;
            cout.width(15); cout << Total_CQ;
            cout.unsetf(ios_base::floatfield);
          }
          else if (equiv_area) { cout.width(15); cout << min(10000.0, max(-10000.0, Total_CLift)); cout.width(15); cout << min(10000.0, max(-10000.0, Total_CDrag)); cout.width(15);
            cout.precision(4);
            cout.setf(ios::scientific, ios::floatfield);
            cout << Total_CNearFieldOF; }
          else if (freesurface) { cout.width(15); cout << Total_CLift; cout.width(15); cout << Total_CFreeSurface; }
          else if (turbo) {
          	cout.setf(ios::scientific, ios::floatfield);
          	switch (config[ZONE_0]->GetKind_TurboPerf(0)) {
          		case BLADE:
          			cout.width(15); cout << KineticEnergyLoss[0]*100.0;
          			cout.width(15); cout << abs((MassFlowIn[0] + MassFlowOut[0])/MassFlowIn[0])*100.0;
          			break;
          		case STAGE: case TURBINE:
          			cout.width(15); cout << TotalStaticEfficiency[0]*100.0;
          			cout.width(15); cout << PressureOut[0]*config[ZONE_0]->GetPressure_Ref();
          			break;
          		default:
          			break;
          	}
          	cout.unsetf(ios_base::floatfield);
          }
          else { cout.width(15); cout << min(10000.0, max(-10000.0, Total_CLift)); cout.width(15); cout << min(10000.0, max(-10000.0, Total_CDrag)); }
          if (aeroelastic) {
            cout.setf(ios::scientific, ios::floatfield);
            cout.width(15); cout << aeroelastic_plunge[0]; //Only output the first marker being monitored to the console.
            cout.width(15); cout << aeroelastic_pitch[0];
            cout.unsetf(ios_base::floatfield);
          }
          cout << endl;
          
          break;
          
        case RANS :
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << direct_coeff << flow_resid << turb_resid;
            if (aeroelastic) ConvHist_file[0] << aeroelastic_coeff;
            if (output_per_surface) ConvHist_file[0] << monitoring_coeff;
            if (output_1d) ConvHist_file[0] << oneD_outputs;
            if (output_massflow and !output_1d) ConvHist_file[0] << massflow_outputs;
            if (direct_diff != NO_DERIVATIVE) ConvHist_file[0] << d_direct_coeff;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed, ios::floatfield);
          
          if (incompressible || freesurface) cout.width(13);
         else  cout.width(14);
         cout << log10(residual_flow[0]);
//          else  cout.width(14),
//                 cout << log10(residual_flow[0]),
//                 cout.width(14);
//          if ( nDim==2 ) cout << log10(residual_flow[3]);
//          if ( nDim==3 ) cout << log10(residual_flow[4]);
          
          switch(nVar_Turb) {
            case 1: cout.width(14); cout << log10(residual_turbulent[0]); break;
            case 2: cout.width(14); cout << log10(residual_turbulent[0]);
              cout.width(15); cout << log10(residual_turbulent[1]); break;
          }
          
          if (transition) { cout.width(14); cout << log10(residual_transition[0]); cout.width(14); cout << log10(residual_transition[1]); }
          
          if (rotating_frame && nDim == 3 && !turbo ) {
            cout.setf(ios::scientific, ios::floatfield);
            cout.width(15); cout << Total_CT; cout.width(15);
            cout << Total_CQ;
            cout.unsetf(ios_base::floatfield);
          }
          else if (equiv_area) { cout.width(15); cout << min(10000.0, max(-10000.0, Total_CLift)); cout.width(15); cout << min(10000.0, max(-10000.0, Total_CDrag)); cout.width(15);
            cout.precision(4);
            cout.setf(ios::scientific, ios::floatfield);
            cout << Total_CNearFieldOF; }
          else if (turbo) {
          cout.setf(ios::scientific, ios::floatfield);
        	switch (config[ZONE_0]->GetKind_TurboPerf(0)) {
        		case BLADE:
        			cout.width(15); cout << KineticEnergyLoss[0]*100.0;
        			cout.width(15); cout << abs((MassFlowIn[0] + MassFlowOut[0])/MassFlowIn[0])*100.0;
        			break;
        		case STAGE: case TURBINE:
        			cout.width(15); cout << TotalStaticEfficiency[0]*100.0;
					cout.width(15); cout << PressureOut[0]*config[ZONE_0]->GetPressure_Ref();
					break;
        		default:
        			break;
        	}
		    cout.unsetf(ios_base::floatfield);
		  }
          else { cout.width(15); cout << min(10000.0, max(-10000.0, Total_CLift)); cout.width(15); cout << min(10000.0, max(-10000.0, Total_CDrag)); }
          
          if (aeroelastic) {
            cout.setf(ios::scientific, ios::floatfield);
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
            cout.setf(ios::fixed, ios::floatfield);
            cout.width(13); cout << log10(residual_flow[0]);
            cout.width(14); cout << log10(residual_levelset[0]);
            cout.width(15); cout << Total_CLift;
            cout.width(14); cout << Total_CFreeSurface;
            
            cout << endl;
          }
          
          break;
          
        case WAVE_EQUATION:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << wave_coeff << wave_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed, ios::floatfield);
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
          cout.setf(ios::fixed, ios::floatfield);
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
          cout.setf(ios::fixed, ios::floatfield);
          cout.width(15); cout << log10(residual_fea[0]);
          cout.width(15); cout << log10(residual_fea[1]);
          if (nDim == 3) { cout.width(15); cout << log10(residual_fea[2]); }
          cout.precision(4);
          cout.setf(ios::scientific, ios::floatfield);
          cout.width(14); cout << Total_CFEA;
          cout << endl;
          break;
          
        case ADJ_EULER :              case ADJ_NAVIER_STOKES :
        case DISC_ADJ_EULER:          case DISC_ADJ_NAVIER_STOKES:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed, ios::floatfield);
          if (compressible) {
            cout.width(15); cout << log10(residual_adjflow[0]);
            cout.width(15); cout << log10(residual_adjflow[nDim+1]);
          }
          if (incompressible || freesurface) {
            cout.width(17); cout << log10(residual_adjflow[0]);
            cout.width(16); cout << log10(residual_adjflow[1]);
          }

          if (disc_adj){
            cout.precision(4);
            cout.setf(ios::scientific, ios::floatfield);
            cout.width(14); cout << Total_Sens_Press;
            cout.width(14); cout << Total_Sens_Mach;
          }else{
            cout.precision(4);
            cout.setf(ios::scientific, ios::floatfield);
            cout.width(14); cout << Total_Sens_Geo;
            if (output_1d){
              cout.width(14); cout << Total_Sens_BPress;
            }
            else{
              cout.width(14); cout << Total_Sens_Mach;
            }
          }
          cout << endl;
          cout.unsetf(ios_base::floatfield);
          
          if (freesurface) {
            if (!DualTime_Iteration) {
              ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << adj_levelset_resid << end;
              ConvHist_file[0].flush();
            }
            
            cout.precision(6);
            cout.setf(ios::fixed, ios::floatfield);
            cout.width(17); cout << log10(residual_adjflow[0]);
            cout.width(16); cout << log10(residual_adjlevelset[0]);
            cout.precision(3);
            cout.setf(ios::scientific, ios::floatfield);
            cout.width(12); cout << Total_Sens_Geo;
            if (output_1d){
              cout.width(12); cout << Total_Sens_BPress;
            }
            else{
              cout.width(12); cout << Total_Sens_Mach;
            }
            cout.unsetf(ios_base::floatfield);
            cout << endl;
          }
          
          break;
          
        case ADJ_RANS : case DISC_ADJ_RANS:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid;
            if (!config[val_iZone]->GetFrozen_Visc())
              ConvHist_file[0] << adj_turb_resid;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed, ios::floatfield);
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
          if (disc_adj){
            cout.precision(4);
            cout.setf(ios::scientific, ios::floatfield);
            cout.width(14); cout << Total_Sens_Press;
            cout.width(14); cout << Total_Sens_Mach;
          }else{
            cout.precision(4);
            cout.setf(ios::scientific, ios::floatfield);
            cout.width(14); cout << Total_Sens_Geo;
            if (output_1d){
              cout.width(14); cout << Total_Sens_BPress;
            }
            else{
              cout.width(14); cout << Total_Sens_Mach;
            }
          }
          cout << endl;
          cout.unsetf(ios_base::floatfield);
          if (freesurface) {
            if (!DualTime_Iteration) {
              ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << adj_levelset_resid;
              ConvHist_file[0] << end;
              ConvHist_file[0].flush();
            }
            
            cout.precision(6);
            cout.setf(ios::fixed, ios::floatfield);
            cout.width(17); cout << log10(residual_adjflow[0]);
            cout.width(16); cout << log10(residual_adjlevelset[0]);
            
            cout.precision(4);
            cout.setf(ios::scientific, ios::floatfield);
            cout.width(12); cout << Total_Sens_Geo;
            if (output_1d){
              cout.width(14); cout << Total_Sens_BPress;
            }
            else{
              cout.width(14); cout << Total_Sens_Mach;
            }
            cout << endl;
            cout.unsetf(ios_base::floatfield);
          }
          
          break;
          
      }
      cout.unsetf(ios::fixed);
      
      delete [] residual_flow;
      delete [] residual_turbulent;
      delete [] residual_transition;
      delete [] residual_levelset;
      delete [] residual_wave;
      delete [] residual_fea;
      delete [] residual_heat;
      
      delete [] residual_adjflow;
      delete [] residual_adjturbulent;
      delete [] residual_adjlevelset;
      
      delete [] Surface_CLift;
      delete [] Surface_CDrag;
      delete [] Surface_CSideForce;
      delete [] Surface_CEff;
      delete [] Surface_CFx;
      delete [] Surface_CFy;
      delete [] Surface_CFz;
      delete [] Surface_CMx;
      delete [] Surface_CMy;
      delete [] Surface_CMz;
      delete [] aeroelastic_pitch;
      delete [] aeroelastic_plunge;

      delete [] TotalStaticEfficiency;
      delete [] TotalTotalEfficiency;
      delete [] KineticEnergyLoss;
      delete [] TotalPressureLoss;
      delete [] MassFlowIn;
      delete [] MassFlowOut;
      delete [] FlowAngleIn;
      delete [] FlowAngleOut;
      delete [] EulerianWork;
      delete []	TotalEnthalpyIn;
      delete [] PressureRatio;
      delete [] PressureOut;
      delete [] EnthalpyOut;
      delete [] MachIn;
      delete [] MachOut;
      delete [] NormalMachIn;
      delete [] NormalMachOut;
      delete [] VelocityOutIs;
    }
  }
}

void COutput::SetCFL_Number(CSolver ****solver_container, CConfig **config, unsigned short val_iZone) {
  
  su2double CFLFactor = 1.0, power = 1.0, CFL = 0.0, CFLMin = 0.0, CFLMax = 0.0, Div = 1.0, Diff = 0.0, MGFactor[100];
  unsigned short iMesh;
  
  unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
  unsigned long ExtIter = config[val_iZone]->GetExtIter();

  RhoRes_New = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_RMS(0);
  switch( config[val_iZone]->GetKind_Solver()){
    case ADJ_EULER : case ADJ_NAVIER_STOKES: case ADJ_RANS:
      RhoRes_New = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(0);
      break;
  }

  if (RhoRes_New < EPS) RhoRes_New = EPS;
  if (RhoRes_Old < EPS) RhoRes_Old = RhoRes_New;
  
  Div = RhoRes_Old/RhoRes_New;
  Diff = RhoRes_New-RhoRes_Old;
  
  /*--- Compute MG factor ---*/
  
  for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
    if (iMesh == MESH_0) MGFactor[iMesh] = 1.0;
    else MGFactor[iMesh] = MGFactor[iMesh-1] * config[val_iZone]->GetCFL(iMesh)/config[val_iZone]->GetCFL(iMesh-1);
  }
  
  if (Div < 1.0) power = config[val_iZone]->GetCFL_AdaptParam(0);
  else power = config[val_iZone]->GetCFL_AdaptParam(1);
  
  /*--- Detect a stall in the residual ---*/
  
  if ((fabs(Diff) <= RhoRes_New*1E-8) && (ExtIter != 0)) { Div = 0.1; power = config[val_iZone]->GetCFL_AdaptParam(1); }
  
  CFLMin = config[val_iZone]->GetCFL_AdaptParam(2);
  CFLMax = config[val_iZone]->GetCFL_AdaptParam(3);
  
  CFLFactor = pow(Div, power);
  
  for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
    CFL = config[val_iZone]->GetCFL(iMesh);
    CFL *= CFLFactor;
    
    if ((iMesh == MESH_0) && (CFL <= CFLMin)) {
      for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
        config[val_iZone]->SetCFL(iMesh, 1.001*CFLMin*MGFactor[iMesh]);
      }
      break;
    }
    if ((iMesh == MESH_0) && (CFL >= CFLMax)) {
      for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++)
        config[val_iZone]->SetCFL(iMesh, 0.999*CFLMax*MGFactor[iMesh]);
      break;
    }
    
    config[val_iZone]->SetCFL(iMesh, CFL);
    
  }
  
  RhoRes_Old = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_RMS(0);
  switch( config[val_iZone]->GetKind_Solver()){
      case ADJ_EULER : case ADJ_NAVIER_STOKES: case ADJ_RANS:
        RhoRes_Old = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(0);
        break;
    }
  
}


void COutput::SetForces_Breakdown(CGeometry ***geometry,
                                  CSolver ****solver_container,
                                  CConfig **config,
                                  CIntegration ***integration,
                                  unsigned short val_iZone) {
  
  char cstr[200];
  unsigned short iMarker_Monitoring;
  ofstream Breakdown_file;
  
  int rank = MASTER_NODE;
  bool compressible       = (config[val_iZone]->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible     = (config[val_iZone]->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface        = (config[val_iZone]->GetKind_Regime() == FREESURFACE);
  bool unsteady           = (config[val_iZone]->GetUnsteady_Simulation() != NO);
  bool viscous            = config[val_iZone]->GetViscous();
  bool grid_movement      = config[val_iZone]->GetGrid_Movement();
  bool gravity            = config[val_iZone]->GetGravityForce();
  bool turbulent          = config[val_iZone]->GetKind_Solver() == RANS;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
  unsigned short nDim = geometry[val_iZone][FinestMesh]->GetnDim();
  bool flow = ((config[val_iZone]->GetKind_Solver() == EULER) || (config[val_iZone]->GetKind_Solver() == NAVIER_STOKES) ||
               (config[val_iZone]->GetKind_Solver() == RANS));
  
  /*--- Output the mean flow solution using only the master node ---*/
  
  if ((rank == MASTER_NODE) && (flow)) {
    
    cout << "Writing the forces breakdown file." << endl;

    /*--- Initialize variables to store information from all domains (direct solution) ---*/
    
    su2double Total_CLift = 0.0, Total_CDrag = 0.0, Total_CSideForce = 0.0, Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0, Total_CEff = 0.0, Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0,
    Inv_CLift = 0.0, Inv_CDrag = 0.0, Inv_CSideForce = 0.0, Inv_CMx = 0.0, Inv_CMy = 0.0, Inv_CMz = 0.0, Inv_CEff = 0.0, Inv_CFx = 0.0, Inv_CFy = 0.0, Inv_CFz = 0.0,
    *Surface_CLift = NULL, *Surface_CDrag = NULL, *Surface_CSideForce = NULL, *Surface_CEff = NULL, *Surface_CFx = NULL, *Surface_CFy = NULL,  *Surface_CFz = NULL, *Surface_CMx = NULL, *Surface_CMy = NULL, *Surface_CMz = NULL,
    *Surface_CLift_Inv = NULL, *Surface_CDrag_Inv = NULL, *Surface_CSideForce_Inv = NULL, *Surface_CEff_Inv = NULL, *Surface_CFx_Inv = NULL, *Surface_CFy_Inv = NULL,  *Surface_CFz_Inv = NULL, *Surface_CMx_Inv = NULL, *Surface_CMy_Inv = NULL, *Surface_CMz_Inv = NULL;
    time_t now = time(0);
    string dt = ctime(&now); dt[24] = '.';

    /*--- Allocate memory for the coefficients being monitored ---*/
    
    Surface_CLift      = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CDrag      = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CSideForce = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CEff       = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFx        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFy        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFz        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMx        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMy        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMz        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    
    Surface_CLift_Inv      = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CDrag_Inv      = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CSideForce_Inv = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CEff_Inv       = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFx_Inv        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFy_Inv        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFz_Inv        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMx_Inv        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMy_Inv        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMz_Inv        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];

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
    
    /*--- Flow inviscid solution coefficients ---*/
    
    Inv_CLift       = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CLift_Inv();
    Inv_CDrag       = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CDrag_Inv();
    Inv_CSideForce  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CSideForce_Inv();
    Inv_CEff        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CEff_Inv();
    Inv_CMx         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CMx_Inv();
    Inv_CMy         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CMy_Inv();
    Inv_CMz         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CMz_Inv();
    Inv_CFx         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CFx_Inv();
    Inv_CFy         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CFy_Inv();
    Inv_CFz         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CFz_Inv();
    
    /*--- Look over the markers being monitored and get the desired values ---*/
    
    for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
      Surface_CLift[iMarker_Monitoring]      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CLift(iMarker_Monitoring);
      Surface_CDrag[iMarker_Monitoring]      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CDrag(iMarker_Monitoring);
      Surface_CSideForce[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CSideForce(iMarker_Monitoring);
      Surface_CEff[iMarker_Monitoring]       = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CEff(iMarker_Monitoring);
      Surface_CFx[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFx(iMarker_Monitoring);
      Surface_CFy[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFy(iMarker_Monitoring);
      Surface_CFz[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFz(iMarker_Monitoring);
      Surface_CMx[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMx(iMarker_Monitoring);
      Surface_CMy[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMy(iMarker_Monitoring);
      Surface_CMz[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMz(iMarker_Monitoring);
      
      Surface_CLift_Inv[iMarker_Monitoring]      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CLift_Inv(iMarker_Monitoring);
      Surface_CDrag_Inv[iMarker_Monitoring]      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CDrag_Inv(iMarker_Monitoring);
      Surface_CSideForce_Inv[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CSideForce_Inv(iMarker_Monitoring);
      Surface_CEff_Inv[iMarker_Monitoring]       = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CEff_Inv(iMarker_Monitoring);
      Surface_CFx_Inv[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFx_Inv(iMarker_Monitoring);
      Surface_CFy_Inv[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFy_Inv(iMarker_Monitoring);
      Surface_CFz_Inv[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFz_Inv(iMarker_Monitoring);
      Surface_CMx_Inv[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMx_Inv(iMarker_Monitoring);
      Surface_CMy_Inv[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMy_Inv(iMarker_Monitoring);
      Surface_CMz_Inv[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMz_Inv(iMarker_Monitoring);
    }
    
    /*--- Write file name with extension ---*/
    
    string filename = config[val_iZone]->GetBreakdown_FileName();
    strcpy (cstr, filename.data());

    Breakdown_file.open(cstr, ios::out);
    
    Breakdown_file << endl <<"-------------------------------------------------------------------------" << endl;
    Breakdown_file <<"|    ___ _   _ ___                                                      |" << endl;
    Breakdown_file <<"|   / __| | | |_  )   Release 4.0.2  \"Cardinal\"                         |" << endl;
    Breakdown_file <<"|   \\__ \\ |_| |/ /                                                      |" << endl;
    Breakdown_file <<"|   |___/\\___//___|   Suite (Computational Fluid Dynamics Code)         |" << endl;
    Breakdown_file << "|                                                                       |" << endl;
    Breakdown_file << "|   Local date and time: " << dt << "                      |" << endl;
    Breakdown_file <<"-------------------------------------------------------------------------" << endl;
    Breakdown_file << "| SU2 Lead Dev.: Dr. Francisco Palacios, Francisco.D.Palacios@boeing.com|" << endl;
    Breakdown_file << "|                Dr. Thomas D. Economon, economon@stanford.edu          |" << endl;
    Breakdown_file <<"-------------------------------------------------------------------------" << endl;
    Breakdown_file << "| SU2 Developers:                                                       |" << endl;
    Breakdown_file << "| - Prof. Juan J. Alonso's group at Stanford University.                |" << endl;
    Breakdown_file << "| - Prof. Piero Colonna's group at Delft University of Technology.      |" << endl;
    Breakdown_file << "| - Prof. Nicolas R. Gauger's group at Kaiserslautern U. of Technology. |" << endl;
    Breakdown_file << "| - Prof. Alberto Guardone's group at Polytechnic University of Milan.  |" << endl;
    Breakdown_file << "| - Prof. Rafael Palacios' group at Imperial College London.            |" << endl;
    Breakdown_file <<"-------------------------------------------------------------------------" << endl;
    Breakdown_file << "| Copyright (C) 2012-2015 SU2, the open-source CFD code.                |" << endl;
    Breakdown_file << "|                                                                       |" << endl;
    Breakdown_file << "| SU2 is free software; you can redistribute it and/or                  |" << endl;
    Breakdown_file << "| modify it under the terms of the GNU Lesser General Public            |" << endl;
    Breakdown_file << "| License as published by the Free Software Foundation; either          |" << endl;
    Breakdown_file << "| version 2.1 of the License, or (at your option) any later version.    |" << endl;
    Breakdown_file << "|                                                                       |" << endl;
    Breakdown_file << "| SU2 is distributed in the hope that it will be useful,                |" << endl;
    Breakdown_file << "| but WITHOUT ANY WARRANTY; without even the implied warranty of        |" << endl;
    Breakdown_file << "| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |" << endl;
    Breakdown_file << "| Lesser General Public License for more details.                       |" << endl;
    Breakdown_file << "|                                                                       |" << endl;
    Breakdown_file << "| You should have received a copy of the GNU Lesser General Public      |" << endl;
    Breakdown_file << "| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   |" << endl;
    Breakdown_file <<"-------------------------------------------------------------------------" << endl;
    
    Breakdown_file.precision(6);
    
    Breakdown_file << endl << endl <<"Problem definition:" << endl << endl;
    
    if (compressible) {
      if (viscous) {
        Breakdown_file << "Viscous flow: Computing pressure using the ideal gas law" << endl;
        Breakdown_file << "based on the free-stream temperature and a density computed" << endl;
        Breakdown_file << "from the Reynolds number." << endl;
      } else {
        Breakdown_file << "Inviscid flow: Computing density based on free-stream" << endl;
        Breakdown_file << "temperature and pressure using the ideal gas law." << endl;
      }
    }
    
    if (grid_movement) Breakdown_file << "Force coefficients computed using MACH_MOTION." << endl;
    else Breakdown_file << "Force coefficients computed using free-stream values." << endl;
    
    if (incompressible || freesurface) {
      Breakdown_file << "Viscous and Inviscid flow: rho_ref, and vel_ref" << endl;
      Breakdown_file << "are based on the free-stream values, p_ref = rho_ref*vel_ref^2." << endl;
      Breakdown_file << "The free-stream value of the pressure is 0." << endl;
      Breakdown_file << "Mach number: "<< config[val_iZone]->GetMach() << ", computed using the Bulk modulus." << endl;
      Breakdown_file << "Angle of attack (deg): "<< config[val_iZone]->GetAoA() << ", computed using the the free-stream velocity." << endl;
      Breakdown_file << "Side slip angle (deg): "<< config[val_iZone]->GetAoS() << ", computed using the the free-stream velocity." << endl;
      if (viscous) Breakdown_file << "Reynolds number: " << config[val_iZone]->GetReynolds() << ", computed using free-stream values."<< endl;
      Breakdown_file << "Only dimensional computation, the grid should be dimensional." << endl;
    }
    
    Breakdown_file <<"-- Input conditions:"<< endl;
    
    if (compressible) {
      switch (config[val_iZone]->GetKind_FluidModel()) {
          
        case STANDARD_AIR:
          Breakdown_file << "Fluid Model: STANDARD_AIR "<< endl;
          Breakdown_file << "Specific gas constant: " << config[val_iZone]->GetGas_Constant();
          if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.m/kg.K." << endl;
          else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.ft/slug.R." << endl;
          Breakdown_file << "Specific gas constant (non-dim): " << config[val_iZone]->GetGas_ConstantND()<< endl;
          Breakdown_file << "Specific Heat Ratio: 1.4000 "<< endl;
          break;
          
        case IDEAL_GAS:
          Breakdown_file << "Fluid Model: IDEAL_GAS "<< endl;
          Breakdown_file << "Specific gas constant: " << config[val_iZone]->GetGas_Constant() << " N.m/kg.K." << endl;
          Breakdown_file << "Specific gas constant (non-dim): " << config[val_iZone]->GetGas_ConstantND()<< endl;
          Breakdown_file << "Specific Heat Ratio: "<< config[val_iZone]->GetGamma() << endl;
          break;
          
        case VW_GAS:
          Breakdown_file << "Fluid Model: Van der Waals "<< endl;
          Breakdown_file << "Specific gas constant: " << config[val_iZone]->GetGas_Constant() << " N.m/kg.K." << endl;
          Breakdown_file << "Specific gas constant (non-dim): " << config[val_iZone]->GetGas_ConstantND()<< endl;
          Breakdown_file << "Specific Heat Ratio: "<< config[val_iZone]->GetGamma() << endl;
          Breakdown_file << "Critical Pressure:   " << config[val_iZone]->GetPressure_Critical()  << " Pa." << endl;
          Breakdown_file << "Critical Temperature:  " << config[val_iZone]->GetTemperature_Critical() << " K." << endl;
          Breakdown_file << "Critical Pressure (non-dim):   " << config[val_iZone]->GetPressure_Critical() /config[val_iZone]->GetPressure_Ref() << endl;
          Breakdown_file << "Critical Temperature (non-dim) :  " << config[val_iZone]->GetTemperature_Critical() /config[val_iZone]->GetTemperature_Ref() << endl;
          break;
          
        case PR_GAS:
          Breakdown_file << "Fluid Model: Peng-Robinson "<< endl;
          Breakdown_file << "Specific gas constant: " << config[val_iZone]->GetGas_Constant() << " N.m/kg.K." << endl;
          Breakdown_file << "Specific gas constant(non-dim): " << config[val_iZone]->GetGas_ConstantND()<< endl;
          Breakdown_file << "Specific Heat Ratio: "<< config[val_iZone]->GetGamma() << endl;
          Breakdown_file << "Critical Pressure:   " << config[val_iZone]->GetPressure_Critical()  << " Pa." << endl;
          Breakdown_file << "Critical Temperature:  " << config[val_iZone]->GetTemperature_Critical() << " K." << endl;
          Breakdown_file << "Critical Pressure (non-dim):   " << config[val_iZone]->GetPressure_Critical() /config[val_iZone]->GetPressure_Ref() << endl;
          Breakdown_file << "Critical Temperature (non-dim) :  " << config[val_iZone]->GetTemperature_Critical() /config[val_iZone]->GetTemperature_Ref() << endl;
          break;
      }
      
      if (viscous) {
        
        switch (config[val_iZone]->GetKind_ViscosityModel()) {
            
          case CONSTANT_VISCOSITY:
            Breakdown_file << "Viscosity Model: CONSTANT_VISCOSITY  "<< endl;
            Breakdown_file << "Laminar Viscosity: " << config[val_iZone]->GetMu_ConstantND()*config[val_iZone]->GetViscosity_Ref();
            if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << endl;
            else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << endl;
            Breakdown_file << "Laminar Viscosity (non-dim): " << config[val_iZone]->GetMu_ConstantND()<< endl;
            break;
            
          case SUTHERLAND:
            Breakdown_file << "Viscosity Model: SUTHERLAND "<< endl;
            Breakdown_file << "Ref. Laminar Viscosity: " << config[val_iZone]->GetMu_RefND()*config[val_iZone]->GetViscosity_Ref();
            if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << endl;
            else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << endl;
            Breakdown_file << "Ref. Temperature: " << config[val_iZone]->GetMu_Temperature_RefND()*config[val_iZone]->GetTemperature_Ref();
            if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " K." << endl;
            else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " R." << endl;
            Breakdown_file << "Sutherland Constant: "<< config[val_iZone]->GetMu_SND()*config[val_iZone]->GetTemperature_Ref();
            if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " K." << endl;
            else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " R." << endl;
            Breakdown_file << "Laminar Viscosity (non-dim): " << config[val_iZone]->GetMu_ConstantND()<< endl;
            Breakdown_file << "Ref. Temperature (non-dim): " << config[val_iZone]->GetMu_Temperature_RefND()<< endl;
            Breakdown_file << "Sutherland constant (non-dim): "<< config[val_iZone]->GetMu_SND()<< endl;
            break;
            
        }
        switch (config[val_iZone]->GetKind_ConductivityModel()) {
            
          case CONSTANT_PRANDTL:
            Breakdown_file << "Conductivity Model: CONSTANT_PRANDTL  "<< endl;
            Breakdown_file << "Prandtl: " << config[val_iZone]->GetPrandtl_Lam()<< endl;
            break;
            
          case CONSTANT_CONDUCTIVITY:
            Breakdown_file << "Conductivity Model: CONSTANT_CONDUCTIVITY "<< endl;
            Breakdown_file << "Molecular Conductivity: " << config[val_iZone]->GetKt_ConstantND()*config[val_iZone]->GetConductivity_Ref()<< " W/m^2.K." << endl;
            Breakdown_file << "Molecular Conductivity (non-dim): " << config[val_iZone]->GetKt_ConstantND()<< endl;
            break;
            
        }
      }
    }
    
    if (incompressible || freesurface) {
      Breakdown_file << "Bulk modulus: " << config[val_iZone]->GetBulk_Modulus();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << endl;
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << endl;
      Breakdown_file << "Artificial compressibility factor: " << config[val_iZone]->GetArtComp_Factor();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << endl;
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << endl;
    }
    
    Breakdown_file << "Free-stream static pressure: " << config[val_iZone]->GetPressure_FreeStream();
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << endl;
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << endl;
    
    Breakdown_file << "Free-stream total pressure: " << config[val_iZone]->GetPressure_FreeStream() * pow( 1.0+config[val_iZone]->GetMach()*config[val_iZone]->GetMach()*0.5*(config[val_iZone]->GetGamma()-1.0), config[val_iZone]->GetGamma()/(config[val_iZone]->GetGamma()-1.0) );
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << endl;
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << endl;
    
    if (compressible) {
      Breakdown_file << "Free-stream temperature: " << config[val_iZone]->GetTemperature_FreeStream();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " K." << endl;
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " R." << endl;
    }
    
    Breakdown_file << "Free-stream density: " << config[val_iZone]->GetDensity_FreeStream();
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " kg/m^3." << endl;
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " slug/ft^3." << endl;
    
    if (nDim == 2) {
      Breakdown_file << "Free-stream velocity: (" << config[val_iZone]->GetVelocity_FreeStream()[0] << ", ";
      Breakdown_file << config[val_iZone]->GetVelocity_FreeStream()[1] << ")";
    }
    if (nDim == 3) {
      Breakdown_file << "Free-stream velocity: (" << config[val_iZone]->GetVelocity_FreeStream()[0] << ", ";
      Breakdown_file << config[val_iZone]->GetVelocity_FreeStream()[1] << ", " << config[val_iZone]->GetVelocity_FreeStream()[2] << ")";
    }
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m/s. ";
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft/s. ";
    
    Breakdown_file << "Magnitude: "	<< config[val_iZone]->GetModVel_FreeStream();
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m/s." << endl;
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft/s." << endl;
    
    if (compressible) {
      Breakdown_file << "Free-stream total energy per unit mass: " << config[val_iZone]->GetEnergy_FreeStream();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m^2/s^2." << endl;
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft^2/s^2." << endl;
    }
    
    if (viscous) {
      Breakdown_file << "Free-stream viscosity: " << config[val_iZone]->GetViscosity_FreeStream();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << endl;
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << endl;
      if (turbulent) {
        Breakdown_file << "Free-stream turb. kinetic energy per unit mass: " << config[val_iZone]->GetTke_FreeStream();
        if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m^2/s^2." << endl;
        else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft^2/s^2." << endl;
        Breakdown_file << "Free-stream specific dissipation: " << config[val_iZone]->GetOmega_FreeStream();
        if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " 1/s." << endl;
        else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " 1/s." << endl;
      }
    }
    
    if (unsteady) { Breakdown_file << "Total time: " << config[val_iZone]->GetTotal_UnstTime() << " s. Time step: " << config[val_iZone]->GetDelta_UnstTime() << " s." << endl; }
    
    /*--- Print out reference values. ---*/
    
    Breakdown_file <<"-- Reference values:"<< endl;
    
    if (compressible) {
      Breakdown_file << "Reference specific gas constant: " << config[val_iZone]->GetGas_Constant_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.m/kg.K." << endl;
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.ft/slug.R." << endl;
    }
    
    Breakdown_file << "Reference pressure: " << config[val_iZone]->GetPressure_Ref();
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << endl;
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << endl;
    
    if (compressible) {
      Breakdown_file << "Reference temperature: " << config[val_iZone]->GetTemperature_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " K." << endl;
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " R." << endl;
    }
    
    Breakdown_file << "Reference density: " << config[val_iZone]->GetDensity_Ref();
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " kg/m^3." << endl;
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " slug/ft^3." << endl;
    
    Breakdown_file << "Reference velocity: " << config[val_iZone]->GetVelocity_Ref();
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m/s." << endl;
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft/s." << endl;
    
    if (compressible) {
      Breakdown_file << "Reference energy per unit mass: " << config[val_iZone]->GetEnergy_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m^2/s^2." << endl;
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft^2/s^2." << endl;
    }
    
    if (incompressible || freesurface) {
      Breakdown_file << "Reference length: " << config[val_iZone]->GetLength_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m." << endl;
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " in." << endl;
    }
    
    if (viscous) {
      Breakdown_file << "Reference viscosity: " << config[val_iZone]->GetViscosity_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << endl;
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << endl;
      Breakdown_file << "Reference conductivity: " << config[val_iZone]->GetConductivity_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " W/m^2.K." << endl;
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf/ft.s.R." << endl;
    }
    
    
    if (unsteady) Breakdown_file << "Reference time: " << config[val_iZone]->GetTime_Ref() <<" s." << endl;
    
    /*--- Print out resulting non-dim values here. ---*/
    
    Breakdown_file << "-- Resulting non-dimensional state:" << endl;
    Breakdown_file << "Mach number (non-dim): " << config[val_iZone]->GetMach() << endl;
    if (viscous) {
      Breakdown_file << "Reynolds number (non-dim): " << config[val_iZone]->GetReynolds() <<". Re length: " << config[val_iZone]->GetLength_Reynolds();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m." << endl;
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft." << endl;
    }
    if (gravity) {
      Breakdown_file << "Froude number (non-dim): " << config[val_iZone]->GetFroude() << endl;
      Breakdown_file << "Lenght of the baseline wave (non-dim): " << 2.0*PI_NUMBER*config[val_iZone]->GetFroude()*config[val_iZone]->GetFroude() << endl;
    }
    
    if (compressible) {
      Breakdown_file << "Specific gas constant (non-dim): " << config[val_iZone]->GetGas_ConstantND() << endl;
      Breakdown_file << "Free-stream temperature (non-dim): " << config[val_iZone]->GetTemperature_FreeStreamND() << endl;
    }
    
    Breakdown_file << "Free-stream pressure (non-dim): " << config[val_iZone]->GetPressure_FreeStreamND() << endl;
    
    Breakdown_file << "Free-stream density (non-dim): " << config[val_iZone]->GetDensity_FreeStreamND() << endl;
    
    if (nDim == 2) {
      Breakdown_file << "Free-stream velocity (non-dim): (" << config[val_iZone]->GetVelocity_FreeStreamND()[0] << ", ";
      Breakdown_file << config[val_iZone]->GetVelocity_FreeStreamND()[1] << "). ";
    } else {
      Breakdown_file << "Free-stream velocity (non-dim): (" << config[val_iZone]->GetVelocity_FreeStreamND()[0] << ", ";
      Breakdown_file << config[val_iZone]->GetVelocity_FreeStreamND()[1] << ", " << config[val_iZone]->GetVelocity_FreeStreamND()[2] << "). ";
    }
    Breakdown_file << "Magnitude: "	 << config[val_iZone]->GetModVel_FreeStreamND() << endl;
    
    if (compressible)
      Breakdown_file << "Free-stream total energy per unit mass (non-dim): " << config[val_iZone]->GetEnergy_FreeStreamND() << endl;
    
    if (viscous) {
      Breakdown_file << "Free-stream viscosity (non-dim): " << config[val_iZone]->GetViscosity_FreeStreamND() << endl;
      if (turbulent) {
        Breakdown_file << "Free-stream turb. kinetic energy (non-dim): " << config[val_iZone]->GetTke_FreeStreamND() << endl;
        Breakdown_file << "Free-stream specific dissipation (non-dim): " << config[val_iZone]->GetOmega_FreeStreamND() << endl;
      }
    }
    
    if (unsteady) {
      Breakdown_file << "Total time (non-dim): " << config[val_iZone]->GetTotal_UnstTimeND() << endl;
      Breakdown_file << "Time step (non-dim): " << config[val_iZone]->GetDelta_UnstTimeND() << endl;
    }
    
    Breakdown_file << endl << endl <<"Forces breakdown:" << endl << endl;
    
    Breakdown_file << "Total CL:    ";
    Breakdown_file.width(11); Breakdown_file << Total_CLift;
    Breakdown_file << " | Pressure Component (";
    Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Inv_CLift*100.0)/(Total_CLift+EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11); Breakdown_file << Inv_CLift;
    Breakdown_file << " | Friction Component (";
    Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Inv_CLift*100.0)/(Total_CLift+EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11); Breakdown_file << Total_CLift-Inv_CLift << endl;
    
    Breakdown_file << "Total CD:    ";
    Breakdown_file.width(11); Breakdown_file << Total_CDrag;
    Breakdown_file << " | Pressure Component (";
    Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Inv_CDrag*100.0)/(Total_CDrag+EPS)) << "%): ";;
    Breakdown_file.width(11); Breakdown_file << Inv_CDrag;
    Breakdown_file << " | Friction Component (";
    Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Inv_CDrag*100.0)/(Total_CDrag+EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11); Breakdown_file << Total_CDrag-Inv_CDrag << endl;
    
    if (nDim == 3) {
      Breakdown_file << "Total CSF:   ";
      Breakdown_file.width(11); Breakdown_file << Total_CSideForce;
      Breakdown_file << " | Pressure Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Inv_CSideForce*100.0)/(Total_CSideForce+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Inv_CSideForce;
      Breakdown_file << " | Friction Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Inv_CSideForce*100.0)/(Total_CSideForce+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Total_CSideForce-Inv_CSideForce << endl;
    }

    Breakdown_file << "Total CL/CD: ";
    Breakdown_file.width(11); Breakdown_file << Total_CEff;
    Breakdown_file << " | Pressure Component (";
    Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Inv_CEff*100.0)/(Total_CEff+EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11); Breakdown_file << Inv_CEff;
    Breakdown_file << " | Friction Component (";
    Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Inv_CEff*100.0)/(Total_CEff+EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11); Breakdown_file << Total_CEff-Inv_CEff << endl;

    if (nDim == 3) {
      Breakdown_file << "Total CMx:   ";
      Breakdown_file.width(11); Breakdown_file << Total_CMx;
      Breakdown_file << " | Pressure Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Inv_CMx*100.0)/(Total_CMx+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Inv_CMx;
      Breakdown_file << " | Friction Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Inv_CMx*100.0)/(Total_CMx+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Total_CMx-Inv_CMx << endl;
      
      Breakdown_file << "Total CMy:   ";
      Breakdown_file.width(11); Breakdown_file << Total_CMy;
      Breakdown_file << " | Pressure Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Inv_CMy*100.0)/(Total_CMy+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Inv_CMy;
      Breakdown_file << " | Friction Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Inv_CMy*100.0)/(Total_CMy+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Total_CMy-Inv_CMy << endl;
    }

    Breakdown_file << "Total CMz:   ";
    Breakdown_file.width(11); Breakdown_file << Total_CMz;
    Breakdown_file << " | Pressure Component (";
    Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Inv_CMz*100.0)/(Total_CMz+EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11); Breakdown_file << Inv_CMz;
    Breakdown_file << " | Friction Component (";
    Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Inv_CMz*100.0)/(Total_CMz+EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11); Breakdown_file << Total_CMz-Inv_CMz << endl;

    Breakdown_file << "Total CFx:   ";
    Breakdown_file.width(11); Breakdown_file << Total_CFx;
    Breakdown_file << " | Pressure Component (";
    Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Inv_CFx*100.0)/(Total_CFx+EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11); Breakdown_file << Inv_CFx;
    Breakdown_file << " | Friction Component (";
    Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Inv_CFx*100.0)/(Total_CFx+EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11); Breakdown_file << Total_CFx-Inv_CFx << endl;

    Breakdown_file << "Total CFy:   ";
    Breakdown_file.width(11); Breakdown_file << Total_CFy;
    Breakdown_file << " | Pressure Component (";
    Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Inv_CFy*100.0)/(Total_CFy+EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11); Breakdown_file << Inv_CFy;
    Breakdown_file << " | Friction Component (";
    Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Inv_CFy*100.0)/(Total_CFy+EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11); Breakdown_file << Total_CFy-Inv_CFy << endl;

    if (nDim == 3) {
      Breakdown_file << "Total CFz:   ";
      Breakdown_file.width(11); Breakdown_file << Total_CFz;
      Breakdown_file << " | Pressure Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Inv_CFz*100.0)/(Total_CFz+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Inv_CFz;
      Breakdown_file << " | Friction Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Inv_CFz*100.0)/(Total_CFz+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Total_CFz-Inv_CFz << endl;
    }
    
    Breakdown_file << endl << endl;

    for (iMarker_Monitoring = 0; iMarker_Monitoring < config[val_iZone]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
      
      Breakdown_file << "Surface name: " << config[val_iZone]->GetMarker_Monitoring(iMarker_Monitoring) << endl << endl;
      
      Breakdown_file << "Total CL    (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CLift[iMarker_Monitoring]*100.0)/(Total_CLift+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CLift[iMarker_Monitoring];
      Breakdown_file << " | Pressure Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CLift_Inv[iMarker_Monitoring]*100.0)/(Surface_CLift[iMarker_Monitoring]+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CLift_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Surface_CLift_Inv[iMarker_Monitoring]*100.0)/(Surface_CLift[iMarker_Monitoring]+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CLift[iMarker_Monitoring]-Surface_CLift_Inv[iMarker_Monitoring] << endl;

      Breakdown_file << "Total CD    (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CDrag[iMarker_Monitoring]*100.0)/(Total_CDrag+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CDrag[iMarker_Monitoring];
      Breakdown_file << " | Pressure Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CDrag_Inv[iMarker_Monitoring]*100.0)/(Surface_CDrag[iMarker_Monitoring]+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CDrag_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Surface_CDrag_Inv[iMarker_Monitoring]*100.0)/(Surface_CDrag[iMarker_Monitoring]+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CDrag[iMarker_Monitoring]-Surface_CDrag_Inv[iMarker_Monitoring] << endl;
      
      if (nDim == 3) {
        Breakdown_file << "Total CSF   (";
        Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CSideForce[iMarker_Monitoring]*100.0)/(Total_CSideForce+EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11); Breakdown_file << Surface_CSideForce[iMarker_Monitoring];
        Breakdown_file << " | Pressure Component (";
        Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CSideForce_Inv[iMarker_Monitoring]*100.0)/(Surface_CSideForce[iMarker_Monitoring]+EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11); Breakdown_file << Surface_CSideForce_Inv[iMarker_Monitoring];
        Breakdown_file << " | Friction Component (";
        Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Surface_CSideForce_Inv[iMarker_Monitoring]*100.0)/(Surface_CSideForce[iMarker_Monitoring]+EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11); Breakdown_file << Surface_CSideForce[iMarker_Monitoring]-Surface_CSideForce_Inv[iMarker_Monitoring] << endl;
      }
      
      Breakdown_file << "Total CL/CD (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CEff[iMarker_Monitoring]*100.0)/(Total_CEff+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CEff[iMarker_Monitoring];
      Breakdown_file << " | Pressure Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CEff_Inv[iMarker_Monitoring]*100.0)/(Surface_CEff[iMarker_Monitoring]+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CEff_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Surface_CEff_Inv[iMarker_Monitoring]*100.0)/(Surface_CEff[iMarker_Monitoring]+EPS));
      Breakdown_file << "%): ";

      Breakdown_file.width(11); Breakdown_file << Surface_CEff[iMarker_Monitoring]-Surface_CEff_Inv[iMarker_Monitoring] << endl;
      
      if (nDim == 3) {
        
        Breakdown_file << "Total CMx   (";
        Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CMx[iMarker_Monitoring]*100.0)/(Total_CMx+EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11); Breakdown_file << Surface_CMx[iMarker_Monitoring];
        Breakdown_file << " | Pressure Component (";
        Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CMx_Inv[iMarker_Monitoring]*100.0)/(Surface_CMx[iMarker_Monitoring]+EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11); Breakdown_file << Surface_CMx_Inv[iMarker_Monitoring];
        Breakdown_file << " | Friction Component (";
        Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Surface_CMx_Inv[iMarker_Monitoring]*100.0)/(Surface_CMx[iMarker_Monitoring]+EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11); Breakdown_file << Surface_CMx[iMarker_Monitoring]-Surface_CMx_Inv[iMarker_Monitoring] << endl;
        
        Breakdown_file << "Total CMy   (";
        Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CMy[iMarker_Monitoring]*100.0)/(Total_CMy+EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11); Breakdown_file << Surface_CMy[iMarker_Monitoring];
        Breakdown_file << " | Pressure Component (";
        Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CMy_Inv[iMarker_Monitoring]*100.0)/(Surface_CMy[iMarker_Monitoring]+EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11); Breakdown_file << Surface_CMy_Inv[iMarker_Monitoring];
        Breakdown_file << " | Friction Component (";
        Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Surface_CMy_Inv[iMarker_Monitoring]*100.0)/(Surface_CMy[iMarker_Monitoring]+EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11); Breakdown_file << Surface_CMy[iMarker_Monitoring]-Surface_CMy_Inv[iMarker_Monitoring] << endl;
      }
      
      Breakdown_file << "Total CMz   (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CMz[iMarker_Monitoring]*100.0)/(Total_CMz+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CMz[iMarker_Monitoring];
      Breakdown_file << " | Pressure Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CMz_Inv[iMarker_Monitoring]*100.0)/(Surface_CMz[iMarker_Monitoring]+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CMz_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Surface_CMz_Inv[iMarker_Monitoring]*100.0)/(Surface_CMz[iMarker_Monitoring]+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CMz[iMarker_Monitoring]-Surface_CMz_Inv[iMarker_Monitoring] << endl;
      
      Breakdown_file << "Total CFx   (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CFx[iMarker_Monitoring]*100.0)/(Total_CFx+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CFx[iMarker_Monitoring];
      Breakdown_file << " | Pressure Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CFx_Inv[iMarker_Monitoring]*100.0)/(Surface_CFx[iMarker_Monitoring]+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CFx_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Surface_CFx_Inv[iMarker_Monitoring]*100.0)/(Surface_CFx[iMarker_Monitoring]+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CFx[iMarker_Monitoring]-Surface_CFx_Inv[iMarker_Monitoring] << endl;
      
      Breakdown_file << "Total CFy   (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CFy[iMarker_Monitoring]*100.0)/(Total_CFy+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CFy[iMarker_Monitoring];
      Breakdown_file << " | Pressure Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CFy_Inv[iMarker_Monitoring]*100.0)/(Surface_CFy[iMarker_Monitoring]+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CFy_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction Component (";
      Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Surface_CFy_Inv[iMarker_Monitoring]*100.0)/(Surface_CFy[iMarker_Monitoring]+EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11); Breakdown_file << Surface_CFy[iMarker_Monitoring]-Surface_CFy_Inv[iMarker_Monitoring] << endl;
      
      if (nDim == 3) {
        Breakdown_file << "Total CFz   (";
        Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CFz[iMarker_Monitoring]*100.0)/(Total_CFz+EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11); Breakdown_file << Surface_CFz[iMarker_Monitoring];
        Breakdown_file << " | Pressure Component (";
        Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int((Surface_CFz_Inv[iMarker_Monitoring]*100.0)/(Surface_CFz[iMarker_Monitoring]+EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11); Breakdown_file << Surface_CFz_Inv[iMarker_Monitoring];
        Breakdown_file << " | Friction Component (";
        Breakdown_file.width(5); Breakdown_file << SU2_TYPE::Int(100.0-(Surface_CFz_Inv[iMarker_Monitoring]*100.0)/(Surface_CFz[iMarker_Monitoring]+EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11); Breakdown_file << Surface_CFz[iMarker_Monitoring]-Surface_CFz_Inv[iMarker_Monitoring] << endl;
      }
      
      Breakdown_file << endl;

    }
    
    delete [] Surface_CLift;
    delete [] Surface_CDrag;
    delete [] Surface_CSideForce;
    delete [] Surface_CEff;
    delete [] Surface_CFx;
    delete [] Surface_CFy;
    delete [] Surface_CFz;
    delete [] Surface_CMx;
    delete [] Surface_CMy;
    delete [] Surface_CMz;
    
    delete [] Surface_CLift_Inv;
    delete [] Surface_CDrag_Inv;
    delete [] Surface_CSideForce_Inv;
    delete [] Surface_CEff_Inv;
    delete [] Surface_CFx_Inv;
    delete [] Surface_CFy_Inv;
    delete [] Surface_CFz_Inv;
    delete [] Surface_CMx_Inv;
    delete [] Surface_CMy_Inv;
    delete [] Surface_CMz_Inv;
    
    Breakdown_file.close();
    
  }
  
}

void COutput::SetResult_Files(CSolver ****solver_container, CGeometry ***geometry, CConfig **config,
                              unsigned long iExtIter, unsigned short val_nZone) {
  
  int rank = MASTER_NODE;
  
#ifdef HAVE_MPI
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  unsigned short iZone;
  
  for (iZone = 0; iZone < val_nZone; iZone++) {
    
    /*--- Flags identifying the types of files to be written. ---*/
    
    bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
    bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();
    
#ifdef HAVE_MPI
    /*--- Do not merge the volume solutions if we are running in parallel.
     Force the use of SU2_SOL to merge the volume sols in this case. ---*/
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size > SINGLE_NODE) {
      Wrt_Vol = false;
      Wrt_Srf = false;
    }
#endif
    
    bool Wrt_Csv = config[iZone]->GetWrt_Csv_Sol();
    
    if (rank == MASTER_NODE) cout << endl << "Writing comma-separated values (CSV) surface files." << endl;

    switch (config[iZone]->GetKind_Solver()) {
        
      case EULER : case NAVIER_STOKES : case RANS :
        
        if (Wrt_Csv) SetSurfaceCSV_Flow(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][FLOW_SOL], iExtIter, iZone);
        break;
        
      case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS : case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
        if (Wrt_Csv) SetSurfaceCSV_Adjoint(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][ADJFLOW_SOL], solver_container[iZone][MESH_0][FLOW_SOL], iExtIter, iZone);
        break;
        
    }
    
    /*--- Get the file output format ---*/
    
    unsigned short FileFormat = config[iZone]->GetOutput_FileFormat();
    
    /*--- Merge the node coordinates and connectivity, if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/
    
    if (Wrt_Vol || Wrt_Srf) {
      if (rank == MASTER_NODE) cout << "Merging connectivities in the Master node." << endl;
      MergeConnectivity(config[iZone], geometry[iZone][MESH_0], iZone);
    }
    
    /*--- Merge coordinates of all grid nodes (excluding ghost points).
     The grid coordinates are always merged and included first in the
     restart files. ---*/
    
    if (rank == MASTER_NODE) cout << "Merging coordinates in the Master node." << endl;
    MergeCoordinates(config[iZone], geometry[iZone][MESH_0]);

    if ((rank == MASTER_NODE) && (Wrt_Vol || Wrt_Srf)) {
      if (FileFormat == TECPLOT_BINARY) {
        if (rank == MASTER_NODE) cout << "Writing Tecplot binary volume and surface mesh files." << endl;
        SetTecplotBinary_DomainMesh(config[iZone], geometry[iZone][MESH_0], iZone);
        SetTecplotBinary_SurfaceMesh(config[iZone], geometry[iZone][MESH_0], iZone);
        if (!wrote_base_file)
          DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
        if (!wrote_surf_file)
          DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], wrote_surf_file);
      }
    }
    
    /*--- Merge the solution data needed for volume solutions and restarts ---*/
    
    if (rank == MASTER_NODE) cout << "Merging solution in the Master node." << endl;
    MergeSolution(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0], iZone);
    
    /*--- Write restart, or Tecplot files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/
    
    if (rank == MASTER_NODE) {
      
      /*--- Write a native restart file ---*/
      
      if (rank == MASTER_NODE) cout << "Writing SU2 native restart file." << endl;
      SetRestart(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0] , iZone);
      
      if (Wrt_Vol) {
        
        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            
            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file volume solution file." << endl;
            SetTecplotASCII(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
            break;
            
          case FIELDVIEW:
            
            /*--- Write a FieldView ASCII file ---*/
            
            if (rank == MASTER_NODE) cout << "Writing FieldView ASCII file volume solution file." << endl;
            SetFieldViewASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone);
            DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
            break;
            
          case TECPLOT_BINARY:
            
            /*--- Write a Tecplot binary solution file ---*/
            
            if (rank == MASTER_NODE) cout << "Writing Tecplot binary volume solution file." << endl;
            SetTecplotBinary_DomainSolution(config[iZone], geometry[iZone][MESH_0], iZone);
            break;
            
          case FIELDVIEW_BINARY:
            
            /*--- Write a FieldView binary file ---*/
            
            if (rank == MASTER_NODE) cout << "Writing FieldView binary file volume solution file." << endl;
            SetFieldViewBinary(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone);
            DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
            break;

          case PARAVIEW:
            
            /*--- Write a Paraview ASCII file ---*/
            
            if (rank == MASTER_NODE) cout << "Writing Paraview ASCII volume solution file." << endl;
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
            
            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII surface solution file." << endl;
            SetTecplotASCII(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0] , iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
            break;
            
          case TECPLOT_BINARY:
            
            /*--- Write a Tecplot binary solution file ---*/
            
            if (rank == MASTER_NODE) cout << "Writing Tecplot binary surface solution file." << endl;
            SetTecplotBinary_SurfaceSolution(config[iZone], geometry[iZone][MESH_0], iZone);
            break;
            
          case PARAVIEW:
            
            /*--- Write a Paraview ASCII file ---*/
            
            if (rank == MASTER_NODE) cout << "Writing Paraview ASCII surface solution file." << endl;
            SetParaview_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
            break;
            
          default:
            break;
        }
        
      }
      
      /*--- Release memory needed for merging the solution data. ---*/
      
      DeallocateCoordinates(config[iZone], geometry[iZone][MESH_0]);
      DeallocateSolution(config[iZone], geometry[iZone][MESH_0]);
      
    }
    
    /*--- Final broadcast (informing other procs that the base output
     file was written). ---*/
    
#ifdef HAVE_MPI
    SU2_MPI::Bcast(&wrote_base_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&wrote_surf_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif
    
  }
}

void COutput::SetBaselineResult_Files(CSolver **solver, CGeometry **geometry, CConfig **config,
                                      unsigned long iExtIter, unsigned short val_nZone) {

  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  unsigned short iZone;
  
  for (iZone = 0; iZone < val_nZone; iZone++) {
    
    /*--- Flags identifying the types of files to be written. ---*/
    
    bool Low_MemoryOutput = config[iZone]->GetLow_MemoryOutput();
    bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
    bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();
    
    /*--- Get the file output format ---*/
    
    unsigned short FileFormat = config[iZone]->GetOutput_FileFormat();
    
    /*--- Merge the node coordinates and connectivity if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/
    
    if ((Wrt_Vol || Wrt_Srf) && (!Low_MemoryOutput)) {
      if (rank == MASTER_NODE) cout << "Merging connectivities in the Master node." << endl;
      MergeConnectivity(config[iZone], geometry[iZone], iZone);
    }
    
    /*--- Merge the solution data needed for volume solutions and restarts ---*/
    
    if ((Wrt_Vol || Wrt_Srf) && (!Low_MemoryOutput)) {
      if (rank == MASTER_NODE) cout << "Merging solution in the Master node." << endl;
      MergeBaselineSolution(config[iZone], geometry[iZone], solver[iZone], iZone);
    }
    
    /*--- Write restart, Tecplot or Paraview files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/
    
    if (!Low_MemoryOutput) {
      
      if (rank == MASTER_NODE) {
        
        if (Wrt_Vol) {
          
          switch (FileFormat) {
              
            case TECPLOT:
              
              /*--- Write a Tecplot ASCII file ---*/
              
              if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (volume grid)." << endl;
              SetTecplotASCII(config[iZone], geometry[iZone], solver, iZone, val_nZone, false);
              DeallocateConnectivity(config[iZone], geometry[iZone], false);
              break;
              
            case FIELDVIEW:
              
              /*--- Write a FieldView ASCII file ---*/
              
              if (rank == MASTER_NODE) cout << "Writing FieldView ASCII file (volume grid)." << endl;
              SetFieldViewASCII(config[iZone], geometry[iZone], iZone, val_nZone);
              DeallocateConnectivity(config[iZone], geometry[iZone], false);
              break;
              
            case TECPLOT_BINARY:
              
              /*--- Write a Tecplot binary solution file ---*/
              
              if (rank == MASTER_NODE) cout << "Writing Tecplot Binary file (volume grid)." << endl;
              SetTecplotBinary_DomainMesh(config[iZone], geometry[iZone], iZone);
              SetTecplotBinary_DomainSolution(config[iZone], geometry[iZone], iZone);
              break;
              
            case FIELDVIEW_BINARY:
              
              /*--- Write a binary binary file ---*/
              
              if (rank == MASTER_NODE) cout << "Writing FieldView ASCII file (volume grid)." << endl;
              SetFieldViewBinary(config[iZone], geometry[iZone], iZone, val_nZone);
              DeallocateConnectivity(config[iZone], geometry[iZone], false);
              break;

            case PARAVIEW:
              
              /*--- Write a Paraview ASCII file ---*/
              
              if (rank == MASTER_NODE) cout << "Writing Paraview ASCII file (volume grid)." << endl;
              SetParaview_ASCII(config[iZone], geometry[iZone], iZone, val_nZone, false);
              DeallocateConnectivity(config[iZone], geometry[iZone], false);
              break;
              
            default:
              break;
          }
          
        }
        
        if (Wrt_Srf) {
                    
          switch (FileFormat) {
              
            case TECPLOT:
              
              /*--- Write a Tecplot ASCII file ---*/
              
              if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (surface grid)." << endl;
              SetTecplotASCII(config[iZone], geometry[iZone], solver, iZone, val_nZone, true);
              DeallocateConnectivity(config[iZone], geometry[iZone], true);
              break;
              
            case TECPLOT_BINARY:
              
              /*--- Write a Tecplot binary solution file ---*/
              
              if (rank == MASTER_NODE) cout << "Writing Tecplot Binary file (surface grid)." << endl;
              SetTecplotBinary_SurfaceMesh(config[iZone], geometry[iZone], iZone);
              SetTecplotBinary_SurfaceSolution(config[iZone], geometry[iZone], iZone);
              break;
              
            case PARAVIEW:
              
              /*--- Write a Paraview ASCII file ---*/
              
              if (rank == MASTER_NODE) cout << "Writing Paraview ASCII file (surface grid)." << endl;
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
      
    }
  
    else {
      
      if (Wrt_Vol) {
        
        if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (volume grid)." << endl;
        char buffer_char[50], out_file[MAX_STRING_SIZE];
        
        string filename;
        if (!config[iZone]->GetAdjoint()) filename = config[iZone]->GetFlow_FileName();
        else filename = config[iZone]->GetAdj_FileName();

        if (size > 1) {
          SPRINTF (buffer_char, "_%d", SU2_TYPE::Int(rank+1));
          filename = filename + buffer_char;
        }
        
        SPRINTF (buffer_char, ".dat");
        strcpy(out_file, filename.c_str()); strcat(out_file, buffer_char);
        SetTecplotASCII_LowMemory(config[iZone], geometry[iZone], solver, out_file, false);
      }
      
      if (Wrt_Srf) {
        
        if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (surface grid)." << endl;
        char buffer_char[50], out_file[MAX_STRING_SIZE];
        
        string filename;
        if (!config[iZone]->GetAdjoint()) filename = config[iZone]->GetSurfFlowCoeff_FileName();
        else filename = config[iZone]->GetSurfAdjCoeff_FileName();

        if (size > 1) {
          SPRINTF (buffer_char, "_%d", SU2_TYPE::Int(rank+1));
          filename = filename + buffer_char;
        }
        
        SPRINTF (buffer_char, ".dat");
        strcpy(out_file, filename.c_str()); strcat(out_file, buffer_char);
        SetTecplotASCII_LowMemory(config[iZone], geometry[iZone], solver, out_file, true);
      }

    }
    
    /*--- Final broadcast (informing other procs that the base output
     file was written). ---*/
    
#ifdef HAVE_MPI
    SU2_MPI::Bcast(&wrote_base_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif
    
  }
}

void COutput::SetMesh_Files(CGeometry **geometry, CConfig **config, unsigned short val_nZone, bool new_file, bool su2_file) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  unsigned short iZone;
  
  for (iZone = 0; iZone < val_nZone; iZone++) {
    
    /*--- Flags identifying the types of files to be written. ---*/
    
    bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol() && config[iZone]->GetVisualize_Deformation();
    bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol() && config[iZone]->GetVisualize_Deformation();;
    
    /*--- Merge the node coordinates and connectivity if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/
    
    if (rank == MASTER_NODE) cout <<"Merging grid connectivity." << endl;
    MergeConnectivity(config[iZone], geometry[iZone], iZone);
    
    /*--- Merge coordinates of all grid nodes (excluding ghost points).
     The grid coordinates are always merged and included first in the
     restart files. ---*/
    
    if (rank == MASTER_NODE) cout <<"Merging grid coordinates." << endl;
    MergeCoordinates(config[iZone], geometry[iZone]);
    
    /*--- Write restart, Tecplot or Paraview files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/
    
    if (rank == MASTER_NODE) {
      
      if (Wrt_Vol) {
        
        if (rank == MASTER_NODE) cout <<"Writing volume mesh file." << endl;
        
        /*--- Write a Tecplot ASCII file ---*/
        if (config[iZone]->GetOutput_FileFormat()==PARAVIEW) SetParaview_MeshASCII(config[iZone], geometry[iZone], iZone,  val_nZone, false,new_file);
        else SetTecplotASCII_Mesh(config[iZone], geometry[iZone], false, new_file);
        
      }
      
      if (Wrt_Srf) {
        
        if (rank == MASTER_NODE) cout <<"Writing surface mesh file." << endl;
        
        /*--- Write a Tecplot ASCII file ---*/
        if (config[iZone]->GetOutput_FileFormat()==PARAVIEW) SetParaview_MeshASCII(config[iZone], geometry[iZone], iZone,  val_nZone, true,new_file);
        else SetTecplotASCII_Mesh(config[iZone], geometry[iZone], true, new_file);
        

      }

      if (rank == MASTER_NODE) cout <<"Writing .su2 file." << endl;
      
      /*--- Write a .su2 ASCII file ---*/
      
      if (su2_file) SetSU2_MeshASCII(config[iZone], geometry[iZone]);
      
      /*--- Deallocate connectivity ---*/

      DeallocateConnectivity(config[iZone], geometry[iZone], true);
      
    }
    
    /*--- Final broadcast (informing other procs that the base output
     file was written). ---*/
    
#ifdef HAVE_MPI
    SU2_MPI::Bcast(&wrote_base_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif
    
  }
}

void COutput::SetMassFlowRate(CSolver *solver_container, CGeometry *geometry, CConfig *config) {
  unsigned short iDim, iMarker_monitor, iMarker;
  unsigned long iVertex, iPoint;
  su2double Vector[3], Total_Mdot=0.0;
  unsigned short nDim = geometry->GetnDim();

  for (iMarker = 0; iMarker< config->GetnMarker_All(); iMarker++) {
    iMarker_monitor = config->GetMarker_All_Monitoring(iMarker);
    if (iMarker_monitor){
      for (iVertex = 0; iVertex < geometry->nVertex[ iMarker ]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iPoint]->GetDomain()) {
          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);

          for (iDim = 0; iDim < nDim; iDim++)
            Total_Mdot -= Vector[iDim]*(solver_container->node[iPoint]->GetSolution(iDim+1));
        }
      }
    }
  }

#ifdef HAVE_MPI
  /*--- Add AllBound information using all the nodes ---*/
  su2double My_Total_Mdot    = Total_Mdot;    Total_Mdot = 0.0;
  SU2_MPI::Allreduce(&My_Total_Mdot, &Total_Mdot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  /*--- Set the output: reusing same variable from OneDimensionalOutput code ---*/
  solver_container->SetOneD_MassFlowRate(Total_Mdot);
}

void COutput::OneDimensionalOutput(CSolver *solver_container, CGeometry *geometry, CConfig *config) {
  
  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Out1D;
  su2double *Normal = NULL, Area = 0.0, UnitNormal[3],
  Tot_Pressure, Mach, Temperature, Pressure = 0.0, Velocity2, Enthalpy, RhoUA, U,// local values at each node (Velocity2 = V^2). U = normal velocity
  AveragePt = 0.0, AverageMach = 0.0, AverageTemperature = 0.0, MassFlowRate = 0.0, // Area Averaged value ( sum / A )
  VelocityRef = 0.0, EnthalpyRef = 0.0, DensityRef = 0.0, PressureRef = 0.0; // Flux conserved values. TemperatureRef follows ideal gas
  su2double TotalArea=0.0;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  su2double Gamma = config->GetGamma();
  unsigned short nDim = geometry->GetnDim();
  
  
  /*--- Loop over the markers ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    Out1D = config->GetMarker_All_Out_1D(iMarker);
    
    /*--- Loop over the vertices to compute the output ---*/
    
    
    if (Out1D == YES) {
      
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        /*--- Find the normal direction ---*/
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          
          /*--- Compute area, and unitary normal ---*/
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
          for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;
          
          if (compressible)                    Pressure = solver_container->node[iPoint]->GetPressure();
          if (incompressible || freesurface)   Pressure = solver_container->node[iPoint]->GetPressureInc();
          
          /*-- Find velocity normal to the marked surface/opening --*/
          
          U = 0.0; RhoUA = 0.0;
          for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
            U += UnitNormal[iDim]*solver_container->node[iPoint]->GetVelocity(iDim);
            RhoUA -=Normal[iDim]*solver_container->node[iPoint]->GetSolution(iDim+1);
          }
          
          Enthalpy = solver_container->node[iPoint]->GetEnthalpy();
          Velocity2 = solver_container->node[iPoint]->GetVelocity2();
          Temperature = solver_container->node[iPoint]->GetTemperature();
          
          Mach = (sqrt(Velocity2))/ solver_container->node[iPoint]->GetSoundSpeed();
          //Stag_Pressure = Pressure*pow((1.0+((Gamma-1.0)/2.0)*pow(Mach, 2.0)),( Gamma/(Gamma-1.0) ) );
          Tot_Pressure = Pressure + 0.5*solver_container->node[iPoint]->GetDensity()*Velocity2;

          AveragePt += Tot_Pressure * Area;
          TotalArea += Area;
          AverageMach += Mach*Area;
          PressureRef += Pressure * Area;
          AverageTemperature += Temperature*Area;
          MassFlowRate += RhoUA; // RhoU is rho * vn * Area
          VelocityRef+=RhoUA*U*U; // rho u A
          EnthalpyRef+=RhoUA*Enthalpy;
          
        }
      }

    }
    
  }
  
#ifdef HAVE_MPI
  
  /*--- Add AllBound information using all the nodes ---*/
  
  su2double My_Area                = TotalArea;          TotalArea = 0.0;
  su2double My_AveragePt           = AveragePt;          AveragePt = 0.0;
  su2double My_AverageMach         = AverageMach;        AverageMach = 0.0;
  su2double My_AverageTemperature  = AverageTemperature; AverageTemperature = 0.0;
  su2double My_MassFlowRate        = MassFlowRate;       MassFlowRate = 0.0;
  su2double My_PressureRef         = PressureRef;        PressureRef = 0.0;
  su2double My_VelocityRef         = VelocityRef;        VelocityRef = 0.0;
  su2double My_EnthalpyRef         = EnthalpyRef;        EnthalpyRef = 0.0;
  su2double My_DensityRef          = DensityRef;         DensityRef = 0.0;
  
  SU2_MPI::Allreduce(&My_Area, &TotalArea, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&My_AveragePt, &AveragePt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&My_AverageMach, &AverageMach, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&My_AverageTemperature, &AverageTemperature, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&My_MassFlowRate, &MassFlowRate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&My_PressureRef, &PressureRef, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&My_VelocityRef, &VelocityRef, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&My_EnthalpyRef , &EnthalpyRef , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&My_DensityRef , &DensityRef , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
#endif
  
  /*--- Set the 1D output ---*/
  /*--- DensityRef depends on the final values of other flux avg variables ---*/
  VelocityRef=sqrt(VelocityRef/MassFlowRate);
  PressureRef=PressureRef/TotalArea;
  EnthalpyRef=EnthalpyRef/MassFlowRate;
  DensityRef =PressureRef*Gamma/(Gamma-1)/(EnthalpyRef-0.5*VelocityRef*VelocityRef);

  /*Area averaged values*/
  solver_container->SetOneD_TotalPress(AveragePt/TotalArea);
  solver_container->SetOneD_Mach(AverageMach/TotalArea);
  solver_container->SetOneD_Temp(AverageTemperature/TotalArea);
  solver_container->SetOneD_MassFlowRate(MassFlowRate);

  /*Flux averaged values*/
  solver_container->SetOneD_FluxAvgPress(PressureRef);
  solver_container->SetOneD_FluxAvgDensity(DensityRef);
  solver_container->SetOneD_FluxAvgVelocity(VelocityRef);
  solver_container->SetOneD_FluxAvgEntalpy(EnthalpyRef);
  
}

void COutput::SetForceSections(CSolver *solver_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {
  
  short iSection, nSection;
  unsigned long iVertex, iPoint;
  su2double *Plane_P0, *Plane_Normal, MinPlane, MaxPlane, *CPressure, MinXCoord, MaxXCoord, Force[3], ForceInviscid[3],
  MomentInviscid[3] = {0.0,0.0,0.0}, MomentDist[3] = {0.0,0.0,0.0}, RefDensity, RefPressure, RefAreaCoeff, *Velocity_Inf, Gas_Constant, Mach2Vel, Mach_Motion, Gamma, RefVel2 = 0.0, factor, NDPressure, *Origin, RefLengthMoment, Alpha, Beta, CDrag_Inv, CLift_Inv, CMy_Inv;
  vector<su2double> Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil, Pressure_Airfoil;
  string Marker_Tag, Slice_Filename, Slice_Ext;
  ofstream Cp_File;
  unsigned short iDim;
  
  bool grid_movement = config->GetGrid_Movement();
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  
  Plane_P0 = new su2double [3];
  Plane_Normal = new su2double [3];
  CPressure = new su2double[geometry->GetnPoint()];
  
  /*--- Compute some reference quantities and necessary values ---*/
  RefDensity = solver_container->GetDensity_Inf();
  RefPressure = solver_container->GetPressure_Inf();
  RefAreaCoeff = config->GetRefAreaCoeff();
  Velocity_Inf = solver_container->GetVelocity_Inf();
  Gamma = config->GetGamma();
  Origin = config->GetRefOriginMoment(0);
  RefLengthMoment  = config->GetRefLengthMoment();
  Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  Beta             = config->GetAoS()*PI_NUMBER/180.0;
  
  if (grid_movement) {
    Gas_Constant = config->GetGas_ConstantND();
    Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (geometry->GetnDim() == 3) {
    
    /*--- Copy the pressure to an auxiliar structure ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      if (compressible) {
        CPressure[iPoint] = (solver_container->node[iPoint]->GetPressure() - RefPressure)*factor*RefAreaCoeff;
      }
      if (incompressible || freesurface) {
        CPressure[iPoint] = (solver_container->node[iPoint]->GetPressureInc() - RefPressure)*factor*RefAreaCoeff;
      }
    }
    
    nSection = config->GetnSections();
    
    for (iSection = 0; iSection < nSection; iSection++) {
      
      /*--- Read the values from the config file ---*/
      
      MinPlane = config->GetSection_Location(0); MaxPlane = config->GetSection_Location(1);
      MinXCoord = -1E6; MaxXCoord = 1E6;
      
      Plane_Normal[0] = 0.0;    Plane_P0[0] = 0.0;
      Plane_Normal[1] = 0.0;    Plane_P0[1] = 0.0;
      Plane_Normal[2] = 0.0;    Plane_P0[2] = 0.0;
      
      Plane_Normal[config->GetAxis_Orientation()] = 1.0;
      Plane_P0[config->GetAxis_Orientation()] = MinPlane + iSection*(MaxPlane - MinPlane)/su2double(nSection-1);
      
      /*--- Compute the airfoil sections (note that we feed in the Cp) ---*/
      
      geometry->ComputeAirfoil_Section(Plane_P0, Plane_Normal,
                                       MinXCoord, MaxXCoord, CPressure,
                                       Xcoord_Airfoil, Ycoord_Airfoil,
                                       Zcoord_Airfoil, Pressure_Airfoil, true,
                                       config);
      
      if ((rank == MASTER_NODE) && (Xcoord_Airfoil.size() == 0)) {
        cout << "Please check the config file, the section "<< iSection+1 <<" has not been detected." << endl;
      }
      
      /*--- Output the pressure on each section (tecplot format) ---*/
      
      if ((rank == MASTER_NODE) && (Xcoord_Airfoil.size() != 0)) {
        
        /*--- Write Cp at each section ---*/
        
        ofstream Cp_File;
        if (iSection == 0) {
          Cp_File.open("cp_sections.dat", ios::out);
          Cp_File << "TITLE = \"Airfoil sections\"" << endl;
          Cp_File << "VARIABLES = \"X\",\"Y\",\"Z\",\"Cp\"" << endl;
        }
        else Cp_File.open("cp_sections.dat", ios::app);
        
        Cp_File << "ZONE T=\"SECTION_"<< (iSection+1) << "\", NODES= "<< Xcoord_Airfoil.size() << ", ELEMENTS= " << Xcoord_Airfoil.size()-1 << ", DATAPACKING= POINT, ZONETYPE= FELINESEG" << endl;
        
        /*--- Coordinates and pressure value ---*/
        
        if (config->GetSystemMeasurements() == SI) {
          for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
            Cp_File << Xcoord_Airfoil[iVertex] <<" "<< Ycoord_Airfoil[iVertex] <<" "<< Zcoord_Airfoil[iVertex] <<" "<< Pressure_Airfoil[iVertex] <<  endl;
          }
        }
        if (config->GetSystemMeasurements() == US) {
          for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
            Cp_File << Xcoord_Airfoil[iVertex]*12.0 <<" "<< Ycoord_Airfoil[iVertex]*12.0 <<" "<< Zcoord_Airfoil[iVertex]*12.0 <<" "<< Pressure_Airfoil[iVertex] <<  endl;
          }
        }
        
        /*--- Basic conectivity ---*/
        
        for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
          Cp_File << iVertex << "\t" << iVertex+1 << "\n";
        }
        
        Cp_File.close();
        
        
        /*--- Compute load distribution ---*/
        
        ForceInviscid[0] = 0.0; ForceInviscid[1] = 0.0; ForceInviscid[2] = 0.0; MomentInviscid[1] = 0.0;
        
        for (iVertex = 0; iVertex < Xcoord_Airfoil.size()-1; iVertex++) {
          
          NDPressure = 0.5*(Pressure_Airfoil[iVertex]+Pressure_Airfoil[iVertex+1]);
          
          Force[0] = -(Zcoord_Airfoil[iVertex+1] - Zcoord_Airfoil[iVertex])*NDPressure;
          Force[1] = 0.0;
          Force[2] = (Xcoord_Airfoil[iVertex+1] - Xcoord_Airfoil[iVertex])*NDPressure;
          
          ForceInviscid[0] += Force[0];
          ForceInviscid[1] += Force[1];
          ForceInviscid[2] += Force[2];
          
          MomentDist[0] = 0.5*(Xcoord_Airfoil[iVertex] + Xcoord_Airfoil[iVertex+1]) - Origin[0];
          MomentDist[1] = 0.5*(Ycoord_Airfoil[iVertex] + Ycoord_Airfoil[iVertex+1]) - Origin[1];
          MomentDist[2] = 0.5*(Zcoord_Airfoil[iVertex] + Zcoord_Airfoil[iVertex+1]) - Origin[3];
          
          MomentInviscid[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLengthMoment;
          
        }
        
        CLift_Inv = fabs( -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha));
        CDrag_Inv = fabs( ForceInviscid[0]*cos(Alpha)*cos(Beta) + ForceInviscid[1]*sin(Beta) + ForceInviscid[2]*sin(Alpha)*cos(Beta));
        CMy_Inv = MomentInviscid[1];
        
        
        /*--- Write load distribution ---*/
        
        ofstream Load_File;
        if (iSection == 0) {
          Load_File.open("load_distribution.dat", ios::out);
          Load_File << "TITLE = \"Load distribution\"" << endl;
          Load_File << "VARIABLES = \"Y\",\"C<sub>L</sub>\",\"C<sub>D</sub>\",\"C<supb>My</sub>\"" << endl;
          Load_File << "ZONE T=\"Wing load distribution\", NODES= "<< nSection << ", ELEMENTS= " << nSection-1 << ", DATAPACKING= POINT, ZONETYPE= FELINESEG" << endl;
        }
        else Load_File.open("load_distribution.dat", ios::app);
        
        /*--- Coordinates and pressure value ---*/
        
        Load_File << Ycoord_Airfoil[0] <<" "<< CLift_Inv <<" "<< CDrag_Inv  <<" "<< CMy_Inv << endl;
        
        /*--- Basic conectivity ---*/
        
        if (iSection == nSection-1) {
          for (iSection = 1; iSection < nSection; iSection++) {
            Load_File << iSection << "\t" << iSection+1 << "\n";
          }
        }
        
        Load_File.close();
        
        
      }
      
    }
    
    
  }
  
  /*--- Delete dynamically allocated memory ---*/
  
  delete [] Plane_P0;
  delete [] Plane_Normal;
  delete [] CPressure;
  
}

void COutput::SetCp_InverseDesign(CSolver *solver_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {
  
  unsigned short iMarker, icommas, Boundary, iDim;
  unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
  su2double XCoord, YCoord, ZCoord, Pressure, PressureCoeff = 0, Cp, CpTarget, *Normal = NULL, Area, PressDiff;
  bool *PointInDomain;
  string text_line, surfCp_filename;
  ifstream Surface_file;
  char buffer[50], cstr[200];
  
  
  nPointLocal = geometry->GetnPoint();
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nPointGlobal = nPointLocal;
#endif
  
  Point2Vertex = new unsigned long[nPointGlobal][2];
  PointInDomain = new bool[nPointGlobal];
  
  for (iPoint = 0; iPoint < nPointGlobal; iPoint ++)
    PointInDomain[iPoint] = false;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary   = config->GetMarker_All_KindBC(iMarker);
    
    if ((Boundary == EULER_WALL             ) ||
        (Boundary == HEAT_FLUX              ) ||
        (Boundary == ISOTHERMAL             ) ||
        (Boundary == NEARFIELD_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        /*--- The Pressure file uses the global numbering ---*/
        
#ifndef HAVE_MPI
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
#else
        iPoint = geometry->node[geometry->vertex[iMarker][iVertex]->GetNode()]->GetGlobalIndex();
#endif
        
        if (geometry->vertex[iMarker][iVertex]->GetNode() < geometry->GetnPointDomain()) {
          Point2Vertex[iPoint][0] = iMarker;
          Point2Vertex[iPoint][1] = iVertex;
          PointInDomain[iPoint] = true;
          solver_container->SetCPressureTarget(iMarker, iVertex, 0.0);
        }
        
      }
    }
  }
  
  /*--- Prepare to read the surface pressure files (CSV) ---*/
  
  surfCp_filename = "TargetCp";
  strcpy (cstr, surfCp_filename.c_str());
  
  /*--- Write file name with extension if unsteady or steady ---*/
  
  if ((config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) ||
      (config->GetUnsteady_Simulation() == TIME_SPECTRAL)) {
    if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.dat", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.dat",  SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.dat",   SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.dat",    SU2_TYPE::Int(iExtIter));
    if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(iExtIter));
  }
  else
    SPRINTF (buffer, ".dat");
  
  strcat (cstr, buffer);
  
  /*--- Read the surface pressure file ---*/
  
  string::size_type position;
  
  Surface_file.open(cstr, ios::in);
  
  if (!(Surface_file.fail())) {
    
    getline(Surface_file, text_line);
    
    while (getline(Surface_file, text_line)) {
      for (icommas = 0; icommas < 50; icommas++) {
        position = text_line.find( ",", 0 );
        if (position!=string::npos) text_line.erase (position,1);
      }
      stringstream  point_line(text_line);
      
      if (geometry->GetnDim() == 2) point_line >> iPoint >> XCoord >> YCoord >> Pressure >> PressureCoeff;
      if (geometry->GetnDim() == 3) point_line >> iPoint >> XCoord >> YCoord >> ZCoord >> Pressure >> PressureCoeff;
      
      if (PointInDomain[iPoint]) {
        
        /*--- Find the vertex for the Point and Marker ---*/
        
        iMarker = Point2Vertex[iPoint][0];
        iVertex = Point2Vertex[iPoint][1];
        
        solver_container->SetCPressureTarget(iMarker, iVertex, PressureCoeff);
        
      }
      
    }
    
    Surface_file.close();
    
  }
  
  /*--- Compute the pressure difference ---*/
  
  PressDiff = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary   = config->GetMarker_All_KindBC(iMarker);
    
    if ((Boundary == EULER_WALL             ) ||
        (Boundary == HEAT_FLUX              ) ||
        (Boundary == ISOTHERMAL             ) ||
        (Boundary == NEARFIELD_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        
        Cp = solver_container->GetCPressure(iMarker, iVertex);
        CpTarget = solver_container->GetCPressureTarget(iMarker, iVertex);
        
        Area = 0.0;
        for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
          Area += Normal[iDim]*Normal[iDim];
        Area = sqrt(Area);
        
        PressDiff += Area * (CpTarget - Cp) * (CpTarget - Cp);
      }
      
    }
  }
  
#ifdef HAVE_MPI
  su2double MyPressDiff = PressDiff;   PressDiff = 0.0;
  SU2_MPI::Allreduce(&MyPressDiff, &PressDiff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /*--- Update the total Cp difference coeffient ---*/
  
  solver_container->SetTotal_CpDiff(PressDiff);
  
  delete[] Point2Vertex;
  
}

void COutput::SetHeat_InverseDesign(CSolver *solver_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {
  
  unsigned short iMarker, icommas, Boundary, iDim;
  unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
  su2double XCoord, YCoord, ZCoord, PressureCoeff, HeatFlux = 0.0, HeatFluxDiff, HeatFluxTarget, *Normal = NULL, Area,
  Pressure, Cf;
  bool *PointInDomain;
  string text_line, surfHeatFlux_filename;
  ifstream Surface_file;
  char buffer[50], cstr[200];
  
  
  nPointLocal = geometry->GetnPoint();
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nPointGlobal = nPointLocal;
#endif
  
  Point2Vertex = new unsigned long[nPointGlobal][2];
  PointInDomain = new bool[nPointGlobal];
  
  for (iPoint = 0; iPoint < nPointGlobal; iPoint ++)
    PointInDomain[iPoint] = false;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary   = config->GetMarker_All_KindBC(iMarker);
    
    if ((Boundary == EULER_WALL             ) ||
        (Boundary == HEAT_FLUX              ) ||
        (Boundary == ISOTHERMAL             ) ||
        (Boundary == NEARFIELD_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        /*--- The Pressure file uses the global numbering ---*/
        
#ifndef HAVE_MPI
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
#else
        iPoint = geometry->node[geometry->vertex[iMarker][iVertex]->GetNode()]->GetGlobalIndex();
#endif
        
        if (geometry->vertex[iMarker][iVertex]->GetNode() < geometry->GetnPointDomain()) {
          Point2Vertex[iPoint][0] = iMarker;
          Point2Vertex[iPoint][1] = iVertex;
          PointInDomain[iPoint] = true;
          solver_container->SetHeatFluxTarget(iMarker, iVertex, 0.0);
        }
      }
    }
  }
  
  /*--- Prepare to read the surface pressure files (CSV) ---*/
  
  surfHeatFlux_filename = "TargetHeatFlux";
  strcpy (cstr, surfHeatFlux_filename.c_str());
  
  /*--- Write file name with extension if unsteady or steady ---*/
  
  if ((config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) ||
      (config->GetUnsteady_Simulation() == TIME_SPECTRAL)) {
    if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.dat", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.dat",  SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.dat",   SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.dat",    SU2_TYPE::Int(iExtIter));
    if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(iExtIter));
  }
  else
    SPRINTF (buffer, ".dat");
  
  strcat (cstr, buffer);
  
  /*--- Read the surface pressure file ---*/
  
  string::size_type position;
  
  Surface_file.open(cstr, ios::in);
  
  if (!(Surface_file.fail())) {
    
    getline(Surface_file, text_line);
    
    while (getline(Surface_file, text_line)) {
      for (icommas = 0; icommas < 50; icommas++) {
        position = text_line.find( ",", 0 );
        if (position!=string::npos) text_line.erase (position,1);
      }
      stringstream  point_line(text_line);
      
      if (geometry->GetnDim() == 2) point_line >> iPoint >> XCoord >> YCoord >> Pressure >> PressureCoeff >> Cf >> HeatFlux;
      if (geometry->GetnDim() == 3) point_line >> iPoint >> XCoord >> YCoord >> ZCoord >> Pressure >> PressureCoeff >> Cf >> HeatFlux;
      
      if (PointInDomain[iPoint]) {
        
        /*--- Find the vertex for the Point and Marker ---*/
        
        iMarker = Point2Vertex[iPoint][0];
        iVertex = Point2Vertex[iPoint][1];
        
        solver_container->SetHeatFluxTarget(iMarker, iVertex, HeatFlux);
        
      }
      
    }
    
    Surface_file.close();
  }
  
  /*--- Compute the pressure difference ---*/
  
  HeatFluxDiff = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary   = config->GetMarker_All_KindBC(iMarker);
    
    if ((Boundary == EULER_WALL             ) ||
        (Boundary == HEAT_FLUX              ) ||
        (Boundary == ISOTHERMAL             ) ||
        (Boundary == NEARFIELD_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        
        HeatFlux = solver_container->GetHeatFlux(iMarker, iVertex);
        HeatFluxTarget = solver_container->GetHeatFluxTarget(iMarker, iVertex);
        
        Area = 0.0;
        for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
          Area += Normal[iDim]*Normal[iDim];
        Area = sqrt(Area);
        
        HeatFluxDiff += Area * (HeatFluxTarget - HeatFlux) * (HeatFluxTarget - HeatFlux);
        
      }
      
    }
  }
  
#ifdef HAVE_MPI
  su2double MyHeatFluxDiff = HeatFluxDiff;   HeatFluxDiff = 0.0;
  SU2_MPI::Allreduce(&MyHeatFluxDiff, &HeatFluxDiff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /*--- Update the total HeatFlux difference coeffient ---*/
  
  solver_container->SetTotal_HeatFluxDiff(HeatFluxDiff);
  
  delete[] Point2Vertex;
  
}

void COutput::SetEquivalentArea(CSolver *solver_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {
  
  ofstream EquivArea_file, FuncGrad_file;
  unsigned short iMarker = 0, iDim;
  short *AzimuthalAngle = NULL;
  su2double Gamma, auxXCoord, auxYCoord, auxZCoord, InverseDesign = 0.0, DeltaX, Coord_i, Coord_j, jp1Coord, *Coord = NULL, MeanFuntion,
  *Face_Normal = NULL, auxArea, auxPress, Mach, Beta, R_Plane, Pressure_Inf,
  ModVelocity_Inf, Velocity_Inf[3], factor, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL,
  *Pressure = NULL, *FaceArea = NULL, *EquivArea = NULL, *TargetArea = NULL, *NearFieldWeight = NULL,
  *Weight = NULL, jFunction, jp1Function;
  unsigned long jVertex, iVertex, iPoint, nVertex_NearField = 0, auxPoint,
  *IdPoint = NULL, *IdDomain = NULL, auxDomain;
  unsigned short iPhiAngle;
  ofstream NearFieldEA_file; ifstream TargetEA_file;
  
  su2double XCoordBegin_OF = config->GetEA_IntLimit(0);
  su2double XCoordEnd_OF = config->GetEA_IntLimit(1);
  
  unsigned short nDim = geometry->GetnDim();
  su2double AoA = -(config->GetAoA()*PI_NUMBER/180.0);
  su2double EAScaleFactor = config->GetEA_ScaleFactor(); // The EA Obj. Func. should be ~ force based Obj. Func.
  
  int rank = MESH_0;
  
  Mach  = config->GetMach();
  Gamma = config->GetGamma();
  Beta = sqrt(Mach*Mach-1.0);
  R_Plane = fabs(config->GetEA_IntLimit(2));
  Pressure_Inf = config->GetPressure_FreeStreamND();
  Velocity_Inf[0] = config->GetVelocity_FreeStreamND()[0];
  Velocity_Inf[1] = config->GetVelocity_FreeStreamND()[1];
  Velocity_Inf[2] = config->GetVelocity_FreeStreamND()[2];
  ModVelocity_Inf = 0;
  for (iDim = 0; iDim < 3; iDim++)
    ModVelocity_Inf += Velocity_Inf[iDim] * Velocity_Inf[iDim];
  
  factor = 4.0*sqrt(2.0*Beta*R_Plane) / (Gamma*Pressure_Inf*Mach*Mach);
  
#ifndef HAVE_MPI
  
  /*--- Compute the total number of points on the near-field ---*/
  
  nVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
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
  
  Xcoord = new su2double[nVertex_NearField];
  Ycoord = new su2double[nVertex_NearField];
  Zcoord = new su2double[nVertex_NearField];
  AzimuthalAngle = new short[nVertex_NearField];
  IdPoint = new unsigned long[nVertex_NearField];
  IdDomain = new unsigned long[nVertex_NearField];
  Pressure = new su2double[nVertex_NearField];
  FaceArea = new su2double[nVertex_NearField];
  EquivArea = new su2double[nVertex_NearField];
  TargetArea = new su2double[nVertex_NearField];
  NearFieldWeight = new su2double[nVertex_NearField];
  Weight = new su2double[nVertex_NearField];
  
  /*--- Copy the boundary information to an array ---*/
  
  nVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
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
            
            su2double YcoordRot = Ycoord[nVertex_NearField];
            su2double ZcoordRot = Xcoord[nVertex_NearField]*sin(AoA) + Zcoord[nVertex_NearField]*cos(AoA);
            
            /*--- Compute the Azimuthal angle (resolution of degress in the Azimuthal angle)---*/
            
            su2double AngleDouble; short AngleInt;
            AngleDouble = fabs(atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER);
            
            /*--- Fix an azimuthal line due to misalignments of the near-field ---*/
            
            su2double FixAzimuthalLine = config->GetFixAzimuthalLine();
            
            if ((AngleDouble >= FixAzimuthalLine - 0.1) && (AngleDouble <= FixAzimuthalLine + 0.1)) AngleDouble = FixAzimuthalLine - 0.1;
            
            AngleInt = SU2_TYPE::Short(floor(AngleDouble + 0.5));
            if (AngleInt >= 0) AzimuthalAngle[nVertex_NearField] = AngleInt;
            else AzimuthalAngle[nVertex_NearField] = 180 + AngleInt;
          }
          
          if (AzimuthalAngle[nVertex_NearField] <= 60) {
            Pressure[nVertex_NearField] = solver_container->node[iPoint]->GetPressure();
            FaceArea[nVertex_NearField] = fabs(Face_Normal[nDim-1]);
            nVertex_NearField ++;
          }
          
        }
      }
  
#else
  
  int nProcessor;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  unsigned long nLocalVertex_NearField = 0, MaxLocalVertex_NearField = 0;
  int iProcessor;
  
  unsigned long *Buffer_Receive_nVertex = NULL;
  if (rank == MASTER_NODE) {
    Buffer_Receive_nVertex = new unsigned long [nProcessor];
  }
  
  /*--- Compute the total number of points of the near-field ghost nodes ---*/
  
  nLocalVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->node[iPoint]->GetCoord();
        
        if (geometry->node[iPoint]->GetDomain())
          if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0))
            nLocalVertex_NearField ++;
      }
  
  unsigned long *Buffer_Send_nVertex = new unsigned long [1];
  Buffer_Send_nVertex[0] = nLocalVertex_NearField;
  
  /*--- Send Near-Field vertex information --*/
  
  SU2_MPI::Allreduce(&nLocalVertex_NearField, &nVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  delete [] Buffer_Send_nVertex;

  su2double *Buffer_Send_Xcoord = new su2double[MaxLocalVertex_NearField];
  su2double *Buffer_Send_Ycoord = new su2double[MaxLocalVertex_NearField];
  su2double *Buffer_Send_Zcoord = new su2double[MaxLocalVertex_NearField];
  unsigned long *Buffer_Send_IdPoint = new unsigned long [MaxLocalVertex_NearField];
  su2double *Buffer_Send_Pressure = new su2double [MaxLocalVertex_NearField];
  su2double *Buffer_Send_FaceArea = new su2double[MaxLocalVertex_NearField];
  
  su2double *Buffer_Receive_Xcoord = NULL;
  su2double *Buffer_Receive_Ycoord = NULL;
  su2double *Buffer_Receive_Zcoord = NULL;
  unsigned long *Buffer_Receive_IdPoint = NULL;
  su2double *Buffer_Receive_Pressure = NULL;
  su2double *Buffer_Receive_FaceArea = NULL;
  
  if (rank == MASTER_NODE) {
    Buffer_Receive_Xcoord = new su2double[nProcessor*MaxLocalVertex_NearField];
    Buffer_Receive_Ycoord = new su2double[nProcessor*MaxLocalVertex_NearField];
    Buffer_Receive_Zcoord = new su2double[nProcessor*MaxLocalVertex_NearField];
    Buffer_Receive_IdPoint = new unsigned long[nProcessor*MaxLocalVertex_NearField];
    Buffer_Receive_Pressure = new su2double[nProcessor*MaxLocalVertex_NearField];
    Buffer_Receive_FaceArea = new su2double[nProcessor*MaxLocalVertex_NearField];
  }
  
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
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
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
            Buffer_Send_Pressure[nLocalVertex_NearField] = solver_container->node[iPoint]->GetPressure();
            Buffer_Send_FaceArea[nLocalVertex_NearField] = fabs(Face_Normal[nDim-1]);
            nLocalVertex_NearField++;
          }
      }
  
  /*--- Send all the information --*/
  
  SU2_MPI::Gather(Buffer_Send_Xcoord, nBuffer_Xcoord, MPI_DOUBLE, Buffer_Receive_Xcoord, nBuffer_Xcoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Ycoord, nBuffer_Ycoord, MPI_DOUBLE, Buffer_Receive_Ycoord, nBuffer_Ycoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Zcoord, nBuffer_Zcoord, MPI_DOUBLE, Buffer_Receive_Zcoord, nBuffer_Zcoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_IdPoint, nBuffer_IdPoint, MPI_UNSIGNED_LONG, Buffer_Receive_IdPoint, nBuffer_IdPoint, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Pressure, nBuffer_Pressure, MPI_DOUBLE, Buffer_Receive_Pressure, nBuffer_Pressure, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_FaceArea, nBuffer_FaceArea, MPI_DOUBLE, Buffer_Receive_FaceArea, nBuffer_FaceArea, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  delete [] Buffer_Send_Xcoord;
  delete [] Buffer_Send_Ycoord;
  delete [] Buffer_Send_Zcoord;
  delete [] Buffer_Send_IdPoint;
  delete [] Buffer_Send_Pressure;
  delete [] Buffer_Send_FaceArea;

  if (rank == MASTER_NODE) {
    
    Xcoord = new su2double[nVertex_NearField];
    Ycoord = new su2double[nVertex_NearField];
    Zcoord = new su2double[nVertex_NearField];
    AzimuthalAngle = new short[nVertex_NearField];
    IdPoint = new unsigned long[nVertex_NearField];
    IdDomain = new unsigned long[nVertex_NearField];
    Pressure = new su2double[nVertex_NearField];
    FaceArea = new su2double[nVertex_NearField];
    EquivArea = new su2double[nVertex_NearField];
    TargetArea = new su2double[nVertex_NearField];
    NearFieldWeight = new su2double[nVertex_NearField];
    Weight = new su2double[nVertex_NearField];
    
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
          
          su2double YcoordRot = Ycoord[nVertex_NearField];
          su2double ZcoordRot = Xcoord[nVertex_NearField]*sin(AoA) + Zcoord[nVertex_NearField]*cos(AoA);
          
          /*--- Compute the Azimuthal angle ---*/
          
          su2double AngleDouble; short AngleInt;
          AngleDouble = fabs(atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER);
          
          /*--- Fix an azimuthal line due to misalignments of the near-field ---*/
          
          su2double FixAzimuthalLine = config->GetFixAzimuthalLine();
          
          if ((AngleDouble >= FixAzimuthalLine - 0.1) && (AngleDouble <= FixAzimuthalLine + 0.1))
            AngleDouble = FixAzimuthalLine - 0.1;
          
          AngleInt = SU2_TYPE::Short(floor(AngleDouble + 0.5));
          
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
    
    delete [] Buffer_Receive_nVertex;
    
    delete [] Buffer_Receive_Xcoord;
    delete [] Buffer_Receive_Ycoord;
    delete [] Buffer_Receive_Zcoord;
    delete [] Buffer_Receive_IdPoint;
    delete [] Buffer_Receive_Pressure;
    delete [] Buffer_Receive_FaceArea;
    
  }
  
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
    
    vector<vector<su2double> > Xcoord_PhiAngle; Xcoord_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > Ycoord_PhiAngle; Ycoord_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > Zcoord_PhiAngle; Zcoord_PhiAngle.resize(PhiAngleList.size());
    vector<vector<unsigned long> > IdPoint_PhiAngle; IdPoint_PhiAngle.resize(PhiAngleList.size());
    vector<vector<unsigned long> > IdDomain_PhiAngle; IdDomain_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > Pressure_PhiAngle; Pressure_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > FaceArea_PhiAngle; FaceArea_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > EquivArea_PhiAngle; EquivArea_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > TargetArea_PhiAngle; TargetArea_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > NearFieldWeight_PhiAngle; NearFieldWeight_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > Weight_PhiAngle; Weight_PhiAngle.resize(PhiAngleList.size());
    
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
    NearFieldEA_file.open("Equivalent_Area.dat", ios::out);
    NearFieldEA_file << "TITLE = \"Equivalent Area evaluation at each azimuthal angle\"" << endl;
    
    if (config->GetSystemMeasurements() == US)
      NearFieldEA_file << "VARIABLES = \"Height (in) at r="<< R_Plane*12.0 << " in. (cyl. coord. system)\"";
    else
      NearFieldEA_file << "VARIABLES = \"Height (m) at r="<< R_Plane << " m. (cylindrical coordinate system)\"";
    
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      if (config->GetSystemMeasurements() == US)
        NearFieldEA_file << ", \"Equivalent Area (ft<sup>2</sup>), <greek>F</greek>= " << PhiAngleList[iPhiAngle] << " deg.\"";
      else
        NearFieldEA_file << ", \"Equivalent Area (m<sup>2</sup>), <greek>F</greek>= " << PhiAngleList[iPhiAngle] << " deg.\"";
    }
    
    NearFieldEA_file << endl;
    for (iVertex = 0; iVertex < EquivArea_PhiAngle[0].size(); iVertex++) {
      
      su2double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
      su2double XcoordRot_init = Xcoord_PhiAngle[0][0]*cos(AoA) - Zcoord_PhiAngle[0][0]*sin(AoA);
      
      if (config->GetSystemMeasurements() == US)
        NearFieldEA_file << scientific << (XcoordRot - XcoordRot_init) * 12.0;
      else
        NearFieldEA_file << scientific << (XcoordRot - XcoordRot_init);
      
      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
        NearFieldEA_file << scientific << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex];
      }
      
      NearFieldEA_file << endl;
      
    }
    NearFieldEA_file.close();
    
    /*--- Read target equivalent area from the configuration file,
     this first implementation requires a complete table (same as the original
     EA table). so... no interpolation. ---*/
    
    vector<vector<su2double> > TargetArea_PhiAngle_Trans;
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
        vector<su2double> row;
        unsigned short iter = 0;
        
        while (is.good()) {
          string token;
          getline(is, token,',');
          
          istringstream js(token);
          
          su2double data;
          js >> data;
          
          /*--- The first element in the table is the coordinate (in or m)---*/
          
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
    
    su2double PhiFactor = 1.0/su2double(PhiAngleList.size());
    
    /*--- Evaluate the objective function ---*/
    
    InverseDesign = 0;
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
        Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;
        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];
        
        su2double Difference = EquivArea_PhiAngle[iPhiAngle][iVertex]-TargetArea_PhiAngle[iPhiAngle][iVertex];
        su2double percentage = fabs(Difference)*100/fabs(TargetArea_PhiAngle[iPhiAngle][iVertex]);
        
        if ((percentage < 0.1) || (Coord_i < XCoordBegin_OF) || (Coord_i > XCoordEnd_OF)) Difference = 0.0;
        
        InverseDesign += EAScaleFactor*PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex]*Difference*Difference;
        
      }
    
    /*--- Evaluate the weight of the nearfield pressure (adjoint input) ---*/
    
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];
        NearFieldWeight_PhiAngle[iPhiAngle][iVertex] = 0.0;
        for (jVertex = iVertex; jVertex < EquivArea_PhiAngle[iPhiAngle].size(); jVertex++) {
          Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex];
          Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;
          
          su2double Difference = EquivArea_PhiAngle[iPhiAngle][jVertex]-TargetArea_PhiAngle[iPhiAngle][jVertex];
          su2double percentage = fabs(Difference)*100/fabs(TargetArea_PhiAngle[iPhiAngle][jVertex]);
          
          if ((percentage < 0.1) || (Coord_j < XCoordBegin_OF) || (Coord_j > XCoordEnd_OF)) Difference = 0.0;
          
          NearFieldWeight_PhiAngle[iPhiAngle][iVertex] += EAScaleFactor*PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex]*2.0*Difference*factor*sqrt(Coord_j-Coord_i);
        }
      }
    
    /*--- Write the Nearfield pressure at each Azimuthal PhiAngle ---*/
    
    EquivArea_file.precision(15);
    EquivArea_file.open("nearfield_flow.dat", ios::out);
    EquivArea_file << "TITLE = \"Equivalent Area evaluation at each azimuthal angle\"" << endl;
    
    if (config->GetSystemMeasurements() == US)
      EquivArea_file << "VARIABLES = \"Height (in) at r="<< R_Plane*12.0 << " in. (cyl. coord. system)\",\"Equivalent Area (ft<sup>2</sup>)\",\"Target Equivalent Area (ft<sup>2</sup>)\",\"Cp\"" << endl;
    else
      EquivArea_file << "VARIABLES = \"Height (m) at r="<< R_Plane << " m. (cylindrical coordinate system)\",\"Equivalent Area (m<sup>2</sup>)\",\"Target Equivalent Area (m<sup>2</sup>)\",\"Cp\"" << endl;
    
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      EquivArea_file << fixed << "ZONE T= \"<greek>F</greek>=" << PhiAngleList[iPhiAngle] << " deg.\"" << endl;
      for (iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++) {
        
        su2double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
        su2double XcoordRot_init = Xcoord_PhiAngle[0][0]*cos(AoA) - Zcoord_PhiAngle[0][0]*sin(AoA);
        
        if (config->GetSystemMeasurements() == US)
          EquivArea_file << scientific << (XcoordRot - XcoordRot_init) * 12.0;
        else
          EquivArea_file << scientific << (XcoordRot - XcoordRot_init);
        
        EquivArea_file << scientific << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex]
        << ", " << TargetArea_PhiAngle[iPhiAngle][iVertex] << ", " << (Pressure_PhiAngle[iPhiAngle][iVertex]-Pressure_Inf)/Pressure_Inf << endl;
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
      su2double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
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
  
#ifndef HAVE_MPI
  
  /*--- Store the value of the NearField coefficient ---*/
  
  solver_container->SetTotal_CEquivArea(InverseDesign);
  
#else
  
  /*--- Send the value of the NearField coefficient to all the processors ---*/
  
  SU2_MPI::Bcast(&InverseDesign, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- Store the value of the NearField coefficient ---*/
  
  solver_container->SetTotal_CEquivArea(InverseDesign);
  
#endif
  
}

