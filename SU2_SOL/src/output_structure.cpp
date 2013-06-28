/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.2
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/output_structure.hpp"


COutput::COutput(void) {

  /*--- Initialize point and connectivity counters to zero. ---*/
  nGlobal_Poin = 0;
  nGlobal_Elem = 0;
  nGlobal_Tria = 0;
  nGlobal_Quad = 0;
  nGlobal_Tetr = 0;
  nGlobal_Hexa = 0;
  nGlobal_Wedg = 0;
  nGlobal_Pyra = 0;
  nGlobal_Line = 0;

  /*--- Initialize CGNS write flag ---*/
  wrote_CGNS_base = false;

  /*--- Initialize Tecplot write flag ---*/
  wrote_Tecplot_base = false;

}

COutput::~COutput(void) { }

void COutput::MergeGeometry(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
  
  /*--- Merge coordinates of all grid nodes (excluding ghost points). ---*/
  
  MergeCoordinates(config, geometry);
  
  /*--- Merge connectivity for each type of element (excluding halos). ---*/
  
  MergeConnectivity(config, geometry, TRIANGLE    );
  MergeConnectivity(config, geometry, RECTANGLE   );
  MergeConnectivity(config, geometry, TETRAHEDRON );
  MergeConnectivity(config, geometry, HEXAHEDRON  );
  MergeConnectivity(config, geometry, WEDGE       );
  MergeConnectivity(config, geometry, PYRAMID     );
  
  /*--- Update total number of elements after merge. ---*/
  
  nGlobal_Elem = nGlobal_Tria + nGlobal_Quad + nGlobal_Tetr +
  nGlobal_Hexa + nGlobal_Pyra + nGlobal_Wedg;
  
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
  
	/*--- Local variables needed for merging the geometry with MPI. ---*/
  
	unsigned long Buffer_Send_nPoin[1], *Buffer_Recv_nPoin = NULL;
	unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
	unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  
  /*--- Each processor sends its local number of nodes to the master. ---*/
  
	nLocalPoint = geometry->GetnPointDomain();
	Buffer_Send_nPoin[0] = nLocalPoint;
	if (rank == MASTER_NODE) Buffer_Recv_nPoin = new unsigned long[nProcessor];
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
    if (geometry->node[iPoint]->GetDomain()) {
      
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

void COutput::MergeConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- Local variables needed on all processors ---*/
  
	unsigned short NODES_PER_ELEMENT;
  
  unsigned long iPoint, iNode, jNode;
	unsigned long iElem = 0, jElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;
  unsigned long *Conn_Elem;
  
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
  Conn_Elem = new unsigned long[nLocalElem*NODES_PER_ELEMENT];
  
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
       halo cells (periodic BC) ---*/
      if (!isHalo) {
        nElem_Total++;
        for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
          Conn_Elem[jNode] = geometry->elem[iElem]->GetNode(iNode);
          
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
  
  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0, pElem = 0;
  unsigned long MaxLocalElem = 0;
  
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
  
  unsigned long *Buffer_Send_Elem = new unsigned long[nBuffer_Scalar];
  unsigned long *Buffer_Recv_Elem = NULL;
  
  int *Buffer_Send_Halo = new int[MaxLocalElem];
  int *Buffer_Recv_Halo = NULL;
  
  /*--- Prepare the receive buffers on the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new unsigned long[nProcessor*nBuffer_Scalar];
    Buffer_Recv_Halo = new int[nProcessor*MaxLocalElem];
    Conn_Elem = new unsigned long[nProcessor*MaxLocalElem*NODES_PER_ELEMENT];
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
        Buffer_Send_Elem[jNode] = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- Check if this is a halo node. If so, flag this element
         as a halo cell. We will use this later to sort and remove
         any duplicates from the connectivity list. ---*/
        if (!geometry->node[iPoint]->GetDomain())
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
  MPI::COMM_WORLD.Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                         Buffer_Recv_Elem, nBuffer_Scalar, MPI::UNSIGNED_LONG,
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
    
    /*--- Loop for flagging duplicate elements so that they are not
     included in the final connectivity list. ---*/
    kElem = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
        
        /*--- Check if this element was originally marked as a halo. ---*/
        if (Buffer_Recv_Halo[kElem+iElem]) {
          
          /*--- Check all other elements flagged as halos on the
           remaining processors for duplicates (start at iProcessor+1). ---*/
          pElem = (iProcessor+1)*MaxLocalElem;
          for (jProcessor = iProcessor+1; jProcessor < nProcessor; jProcessor++) {
            for (jElem = 0; jElem < Buffer_Recv_nElem[jProcessor]; jElem++) {
              
              /*--- Check if this element was originally marked as a halo. ---*/
              if (Buffer_Recv_Halo[pElem+jElem]) {
                
                /*--- Check for a duplicate by comparing the index of each
                 node in the element. ---*/
                bool isDuplicate = true;
                for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
                  if (Buffer_Recv_Elem[kElem*NODES_PER_ELEMENT+iElem*NODES_PER_ELEMENT+iNode] !=
                      Buffer_Recv_Elem[pElem*NODES_PER_ELEMENT+jElem*NODES_PER_ELEMENT+iNode])
                    isDuplicate = false;
                }
                
                /*--- If we have found a duplicate element, set both the
                 original flag and "write" state booleans to false. In this
                 way, this element will not be found as we continue searching
                 and it will not be written to the connectivity list. ---*/
                if (isDuplicate) {
                  Buffer_Recv_Halo[pElem+jElem] = false;
                  Write_Elem[pElem+jElem] = false;
                }
              }
            }
            pElem = (jProcessor+1)*MaxLocalElem;
          }
        }
      }
      kElem = (iProcessor+1)*MaxLocalElem;
    }
    
    /*--- Store the unique connectivity list for this element type. ---*/
    
    jNode = 0; kNode = 0; jElem = 0; nElem_Total = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
        
        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {
          
          /*--- Increment total count for this element type ---*/
          nElem_Total++;
          
          /*--- Get global index, then loop over each variable and store ---*/
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode];
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

void COutput::MergeBaselineSolution(CConfig *config, CGeometry *geometry, CSolution *solution, unsigned short val_iZone) {
  
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
   for nodes in the domain (ignore periodic/sliding halo nodes). ---*/
  jPoint = 0;
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		if (geometry->node[iPoint]->GetDomain()) {
			/*--- Solution (first, and second system of equations) ---*/
			unsigned short jVar = 0;
			for (iVar = 0; iVar < nVar_Total; iVar++) {
				Data[jVar][jPoint] = solution->node[iPoint]->GetSolution(iVar);
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
    
  /*--- Each processor sends its local number of nodes to the master. ---*/
	//!
	//! TO DO: MPI I/O for writing the solution files.
	//!
	nLocalPoint = geometry->GetnPointDomain();
	Buffer_Send_nPoint[0] = nLocalPoint;
	if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Gather(&Buffer_Send_nPoint, 1, MPI::UNSIGNED_LONG,
                         Buffer_Recv_nPoint, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
	nBuffer_Scalar = MaxLocalPoint;

	//!
	//! TO DO: Here is where the code can be extended to an arbitrary number
	//! of variables specified by the user (name & type), and more
	//! logic needs to be done.
	//!
  
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
			if (geometry->node[iPoint]->GetDomain()) {
        
				/*--- Get this variable into the temporary send buffer. ---*/
				Buffer_Send_Var[jPoint] = solution->node[iPoint]->GetSolution(iVar);
        
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

void COutput::SetBaselineResult_Files(CSolution **solution, CGeometry **geometry, CConfig **config,
                              unsigned long iExtIter, unsigned short val_nZone) {
    
    int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
    
	unsigned short iZone;
    
	for (iZone = 0; iZone < val_nZone; iZone++) {

        /*--- Construction/Testing of new I/O routines. The conditional guards will
         come off soon and these will be the default options. ---*/
        
        /*--- Merge the node coordinates and connectivity. ---*/
        if (config[iZone]->GetWrt_Sol_Tec_ASCII()  ||
            config[iZone]->GetWrt_Sol_Tec_Binary() ||
            config[iZone]->GetWrt_Sol_CGNS())
            MergeGeometry(config[iZone], geometry[iZone], iZone);
        
        if (config[iZone]->GetWrt_Sol_Tec_ASCII()  ||
            config[iZone]->GetWrt_Sol_Tec_Binary() ||
            config[iZone]->GetWrt_Sol_CGNS() ||
            config[iZone]->GetWrt_Restart())
            MergeBaselineSolution(config[iZone], geometry[iZone], solution[iZone], iZone);
        
        /*--- Write restart, CGNS, or Tecplot files using the merged data.
         This data lives only on the master, and these routines are currently
         executed by the master proc alone (as if in serial). ---*/
    
        if (rank == MASTER_NODE) {

            /*--- Write a Tecplot ASCII file ---*/
            if (config[iZone]->GetWrt_Sol_Tec_ASCII())
                WriteTecplotASCII(config[iZone], geometry[iZone], iZone, val_nZone);
            
            /*--- Write a Tecplot binary file ---*/
            if (config[iZone]->GetWrt_Sol_Tec_Binary()) WriteTecplotBinary(config[iZone], geometry[iZone], iZone);
            
            /*--- Write a CGNS file ---*/
            if (config[iZone]->GetWrt_Sol_CGNS()) WriteCGNS(config[iZone], geometry[iZone], iZone);
            
            /*--- Clean up memory after writing all output ---*/
            CleanUp(config[iZone], geometry[iZone]);
            
        }
        
    }
}
