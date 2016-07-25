/*!
 * \file definition_structure.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 4.2.0 "Cardinal"
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
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/definition_structure.hpp"

void Driver_Preprocessing(CDriver **driver,
                          CIteration **iteration_container,
                          CSolver ****solver_container,
                          CGeometry ***geometry_container,
                          CIntegration ***integration_container,
                          CNumerics *****numerics_container,
                          CInterpolator ***interpolator_container,
                          CTransfer ***transfer_container,
                          CConfig **config_container,
                          unsigned short val_nZone,
                          unsigned short val_nDim) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- fsi implementations will use, as of now, BGS implentation. More to come. ---*/
  bool fsi = config_container[ZONE_0]->GetFSI_Simulation();
  
  if (val_nZone == SINGLE_ZONE && !config_container[ZONE_0]->GetDiscrete_Adjoint()) {
    
    /*--- Single zone problem: instantiate the single zone driver class. ---*/
    if (rank == MASTER_NODE) cout << "Instantiating a single zone driver for the problem. " << endl;
    
    *driver = new CSingleZoneDriver(iteration_container, solver_container, geometry_container,
                                    integration_container, numerics_container, interpolator_container,
                                    transfer_container, config_container, val_nZone, val_nDim);
    
  } else if (config_container[ZONE_0]->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    
    /*--- Use the spectral method driver. ---*/
    
    if (rank == MASTER_NODE) cout << "Instantiating a spectral method driver for the problem. " << endl;
    *driver = new CSpectralDriver(iteration_container, solver_container, geometry_container,
            					  integration_container, numerics_container, interpolator_container,
            					  transfer_container, config_container, val_nZone, val_nDim);

  } else if ((val_nZone == 2) && fsi) {

	    /*--- FSI problem: instantiate the FSI driver class. ---*/

	 if (rank == MASTER_NODE) cout << "Instantiating a Fluid-Structure Interaction driver for the problem. " << endl;
	 *driver = new CFSIDriver(iteration_container, solver_container, geometry_container,
                      integration_container, numerics_container, interpolator_container,
                            transfer_container, config_container, val_nZone, val_nDim);

  } else if (config_container[ZONE_0]->GetDiscrete_Adjoint()){

    if (rank == MASTER_NODE) cout << "Instantiating a multi-zone discrete adjoint driver for the problem. " << endl;
    *driver = new CDiscAdjMultiZoneDriver(iteration_container, solver_container, geometry_container,
	            			  integration_container, numerics_container, interpolator_container,
	                          transfer_container, config_container, val_nZone, val_nDim);

  } else {
    
    /*--- Multi-zone problem: instantiate the multi-zone driver class by default
     or a specialized driver class for a particular multi-physics problem. ---*/
    
    if (rank == MASTER_NODE) cout << "Instantiating a multi-zone driver for the problem. " << endl;
    *driver = new CMultiZoneDriver(iteration_container, solver_container, geometry_container,
            					   integration_container, numerics_container, interpolator_container,
                                   transfer_container, config_container, val_nZone, val_nDim);
    

    
  }

      /*--- Future multi-zone drivers instatiated here. ---*/
}


void Geometrical_Preprocessing(CGeometry ***geometry, CConfig **config, unsigned short val_nZone) {
  
  unsigned short iMGlevel, iZone;
  unsigned long iPoint;
  int rank = MASTER_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  for (iZone = 0; iZone < val_nZone; iZone++) {
    
    /*--- Compute elements surrounding points, points surrounding points ---*/
    
    if (rank == MASTER_NODE) cout << "Setting point connectivity." << endl;
    geometry[iZone][MESH_0]->SetPoint_Connectivity();
    
    /*--- Renumbering points using Reverse Cuthill McKee ordering ---*/
    
    if (rank == MASTER_NODE) cout << "Renumbering points (Reverse Cuthill McKee Ordering)." << endl;
    geometry[iZone][MESH_0]->SetRCM_Ordering(config[iZone]);
    
    /*--- recompute elements surrounding points, points surrounding points ---*/
    
    if (rank == MASTER_NODE) cout << "Recomputing point connectivity." << endl;
    geometry[iZone][MESH_0]->SetPoint_Connectivity();
    
    /*--- Compute elements surrounding elements ---*/
    
    if (rank == MASTER_NODE) cout << "Setting element connectivity." << endl;
    geometry[iZone][MESH_0]->SetElement_Connectivity();
    
    /*--- Check the orientation before computing geometrical quantities ---*/
    
    if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." << endl;
    geometry[iZone][MESH_0]->SetBoundVolume();
    geometry[iZone][MESH_0]->Check_IntElem_Orientation(config[iZone]);
    geometry[iZone][MESH_0]->Check_BoundElem_Orientation(config[iZone]);
    
    /*--- Create the edge structure ---*/
    
    if (rank == MASTER_NODE) cout << "Identifying edges and vertices." << endl;
    geometry[iZone][MESH_0]->SetEdges();
    geometry[iZone][MESH_0]->SetVertex(config[iZone]);
    
    /*--- Compute cell center of gravity ---*/
    
    if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl;
    geometry[iZone][MESH_0]->SetCoord_CG();
    
    /*--- Create the control volume structures ---*/
    
    if (rank == MASTER_NODE) cout << "Setting the control volume structure." << endl;
    geometry[iZone][MESH_0]->SetControlVolume(config[iZone], ALLOCATE);
    geometry[iZone][MESH_0]->SetBoundControlVolume(config[iZone], ALLOCATE);
    
    /*--- Visualize a dual control volume if requested ---*/
    
    if ((config[iZone]->GetVisualize_CV() >= 0) &&
        (config[iZone]->GetVisualize_CV() < (long)geometry[iZone][MESH_0]->GetnPointDomain()))
      geometry[iZone][MESH_0]->VisualizeControlVolume(config[iZone], UPDATE);
    
    /*--- Identify closest normal neighbor ---*/
    
    if (rank == MASTER_NODE) cout << "Searching for the closest normal neighbors to the surfaces." << endl;
    geometry[iZone][MESH_0]->FindNormal_Neighbor(config[iZone]);
    
    /*--- Compute the surface curvature ---*/
    
    if (rank == MASTER_NODE) cout << "Compute the surface curvature." << endl;
    geometry[iZone][MESH_0]->ComputeSurf_Curvature(config[iZone]);
    
    if ((config[iZone]->GetnMGLevels() != 0) && (rank == MASTER_NODE))
      cout << "Setting the multigrid structure." << endl;
    
    /*--- Create turbovertex structure ---*/
    if (config[iZone]->GetBoolTurbomachinery()){
    	geometry[iZone][MESH_0]->SetTurboVertex(config[iZone], iZone, INFLOW, true);
			geometry[iZone][MESH_0]->SetTurboVertex(config[iZone], iZone, OUTFLOW, true);
    }
	}
  
  /*--- Loop over all the new grid ---*/
  
  for (iMGlevel = 1; iMGlevel <= config[ZONE_0]->GetnMGLevels(); iMGlevel++) {
    
    /*--- Loop over all zones at each grid level. ---*/
    
    for (iZone = 0; iZone < val_nZone; iZone++) {
      
      /*--- Create main agglomeration structure ---*/
      
      geometry[iZone][iMGlevel] = new CMultiGridGeometry(geometry, config, iMGlevel, iZone);
      
      /*--- Compute points surrounding points. ---*/
      
      geometry[iZone][iMGlevel]->SetPoint_Connectivity(geometry[iZone][iMGlevel-1]);
      
      /*--- Create the edge structure ---*/
      
      geometry[iZone][iMGlevel]->SetEdges();
      geometry[iZone][iMGlevel]->SetVertex(geometry[iZone][iMGlevel-1], config[iZone]);
      
      /*--- Create the control volume structures ---*/
      
      geometry[iZone][iMGlevel]->SetControlVolume(config[iZone], geometry[iZone][iMGlevel-1], ALLOCATE);
      geometry[iZone][iMGlevel]->SetBoundControlVolume(config[iZone], geometry[iZone][iMGlevel-1], ALLOCATE);
      geometry[iZone][iMGlevel]->SetCoord(geometry[iZone][iMGlevel-1]);
      
      /*--- Find closest neighbor to a surface point ---*/
      
      geometry[iZone][iMGlevel]->FindNormal_Neighbor(config[iZone]);
      

    }
    
  }
  
  /*--- For unsteady simulations, initialize the grid volumes
   and coordinates for previous solutions. Loop over all zones/grids ---*/
  
  for (iZone = 0; iZone < val_nZone; iZone++) {
    if (config[iZone]->GetUnsteady_Simulation() && config[iZone]->GetGrid_Movement()) {
      for (iMGlevel = 0; iMGlevel <= config[iZone]->GetnMGLevels(); iMGlevel++) {
        for (iPoint = 0; iPoint < geometry[iZone][iMGlevel]->GetnPoint(); iPoint++) {
          
          /*--- Update cell volume ---*/
          
          geometry[iZone][iMGlevel]->node[iPoint]->SetVolume_n();
          geometry[iZone][iMGlevel]->node[iPoint]->SetVolume_nM1();
          
          /*--- Update point coordinates ---*/
          geometry[iZone][iMGlevel]->node[iPoint]->SetCoord_n();
          geometry[iZone][iMGlevel]->node[iPoint]->SetCoord_n1();
          
        }
      }
    }
  }
  
}

void Partition_Analysis(CGeometry *geometry, CConfig *config) {
  
  /*--- This routine does a quick and dirty output of the total
   vertices, ghost vertices, total elements, ghost elements, etc.,
   so that we can analyze the partition quality. ---*/
  
  unsigned short nMarker = config->GetnMarker_All();
  unsigned short iMarker, iNodes, MarkerS, MarkerR;
  unsigned long iElem, iPoint, nVertexS, nVertexR;
  unsigned long nNeighbors = 0, nSendTotal = 0, nRecvTotal = 0;
  unsigned long nPointTotal=0, nPointGhost=0, nElemTotal=0;
  unsigned long nElemHalo=0, nEdge=0, nElemBound=0;
  int iRank;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  nPointTotal = geometry->GetnPoint();
  nPointGhost = geometry->GetnPoint() - geometry->GetnPointDomain();
  nElemTotal  = geometry->GetnElem();
  nEdge       = geometry->GetnEdge();
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    nElemBound  += geometry->GetnElem_Bound(iMarker);
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      MarkerS = iMarker;  MarkerR = iMarker+1;
      nVertexS = geometry->nVertex[MarkerS];
      nVertexR = geometry->nVertex[MarkerR];
      nNeighbors++;
      nSendTotal += nVertexS;
      nRecvTotal += nVertexR;
    }
  }
  
  bool *isHalo = new bool[geometry->GetnElem()];
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    isHalo[iElem] = false;
    for (iNodes = 0; iNodes < geometry->elem[iElem]->GetnNodes(); iNodes++) {
      iPoint = geometry->elem[iElem]->GetNode(iNodes);
      if (!geometry->node[iPoint]->GetDomain()) isHalo[iElem] = true;
    }
  }
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if (isHalo[iElem]) nElemHalo++;
  }
  
  unsigned long *row_ptr = NULL, nnz;
  unsigned short *nNeigh = NULL;
  vector<unsigned long>::iterator it;
  vector<unsigned long> vneighs;
  
  /*--- Don't delete *row_ptr, *col_ind because they are
   asigned to the Jacobian structure. ---*/
  
  /*--- Compute the number of neighbors ---*/
  
  nNeigh = new unsigned short [geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    // +1 -> to include diagonal element
    nNeigh[iPoint] = (geometry->node[iPoint]->GetnPoint()+1);
  }
  
  /*--- Create row_ptr structure, using the number of neighbors ---*/
  
  row_ptr = new unsigned long [geometry->GetnPoint()+1];
  row_ptr[0] = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    row_ptr[iPoint+1] = row_ptr[iPoint] + nNeigh[iPoint];
  nnz = row_ptr[geometry->GetnPoint()];
  
  delete [] row_ptr;
  delete [] nNeigh;
  
  /*--- Now put this info into a CSV file for processing ---*/
  
  char cstr[200];
  ofstream Profile_File;
  strcpy (cstr, "partitioning.csv");
  Profile_File.precision(15);
  
  if (rank == MASTER_NODE) {
    /*--- Prepare and open the file ---*/
    Profile_File.open(cstr, ios::out);
    /*--- Create the CSV header ---*/
    Profile_File << "\"Rank\", \"nNeighbors\", \"nPointTotal\", \"nEdge\", \"nPointGhost\", \"nSendTotal\", \"nRecvTotal\", \"nElemTotal\", \"nElemBoundary\", \"nElemHalo\", \"nnz\"" << endl;
    Profile_File.close();
  }
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- Loop through the map and write the results to the file ---*/
  
  for (iRank = 0; iRank < size; iRank++) {
    if (rank == iRank) {
      Profile_File.open(cstr, ios::out | ios::app);
      Profile_File << rank << ", " << nNeighbors << ", " << nPointTotal << ", " << nEdge << "," << nPointGhost << ", " << nSendTotal << ", " << nRecvTotal << ", " << nElemTotal << "," << nElemBound << ", " << nElemHalo << ", " << nnz << endl;
      Profile_File.close();
    }
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  
  delete [] isHalo;
  
}
