/*!
 * \file SU2_SMC.cpp
 * \brief Main file of Sliding Mesh Code (SU2_SMC).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
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

#include "../include/SU2_SMC.hpp"
using namespace std;

int main(int argc, char *argv[]) {
	
	char file_name[100];
	string MeshFile;
	unsigned short nDim, iZone, nZone = 1;
  
  /*--- Pointer to different structures needed by SU2_SMC ---*/
	CGeometry **geometry_container = NULL;
	CConfig **config_container = NULL;
	
	/*--- Definition of the containers per zone ---*/
	config_container = new CConfig*[MAX_ZONES];
	geometry_container = new CGeometry*[MAX_ZONES];
	
	/*--- Check the number of zones in the mesh file ---*/
	CConfig *config = NULL;
	if (argc == 2) config = new CConfig(argv[1]);
	else { strcpy (file_name, "default.cfg"); config = new CConfig(file_name); }
	nZone = GetnZone(config->GetMesh_FileName(), config->GetMesh_FileFormat(), config);
	nDim  = GetnDim(config->GetMesh_FileName(), config->GetMesh_FileFormat());
    
  for (iZone = 0; iZone < nZone; iZone++) {
		
		/*--- Definition of the configuration class per zones ---*/
		if (argc == 2) config_container[iZone] = new CConfig(argv[1], SU2_SMC, iZone, nZone, VERB_HIGH);
		else {
      strcpy (file_name, "default.cfg");
      config_container[iZone] = new CConfig(file_name, SU2_SMC, iZone, nZone, VERB_HIGH);
    }
  
		/*--- Definition of the geometry class and open the mesh file ---*/
		geometry_container[iZone] = new CPhysicalGeometry(config_container[iZone], config_container[iZone]->GetMesh_FileName(),
                                                      config_container[iZone]->GetMesh_FileFormat(), iZone+1, nZone);
    
    /*--- Perform the geometric preprocessing for each zone ---*/
    Geometric_Definition(geometry_container[iZone], config_container[iZone]);
    
  }
  
  /*--- Begin halo node construction for sliding mesh interfaces. ---*/
  cout << endl <<"--------------------- Sliding Interface Construction --------------------" << endl;
  
   /*--- Match the markers of any sliding interfaces point-by-point. ---*/
  MatchSliding_Interfaces(geometry_container, config_container, nZone);
  
  //!
  //! TO DO: Need to consider overlapping points in the matching phase bc
  //! we only have room to store one 'donor' point at each node currently.
  //! Could perhaps do everything in here one zone at a time...
  //!
  
  //! I have removed lines with this because of the new config structure in each zone:
  //if (config_container[iZone]->GetSlideBound_Zone(config_container[iZone]->GetMarker_All_Tag(iMarker)) == iZone)

  
  /*--- This section finds and replicates all nodes and elements that must be copied. ---*/
  unsigned long iPoint, jPoint, kPoint, pPoint, iVertex, jVertex, iPointSliding, newNodes[4];
  unsigned long iElem, jElem, kElem;
  unsigned short iMarker, jMarker, kMarker, iSliding, iNeighbor, jNeighbor, iNode;
  unsigned short jZone, VTK_Type;
  
  vector<unsigned long>::iterator IterElem[MAX_NUMBER_SLIDING];
  vector<unsigned long>::iterator IterPoint[MAX_NUMBER_SLIDING][2];
  vector<unsigned long>::iterator IterNewElem[MAX_NUMBER_MARKER];
  
  vector<unsigned long> SlidingPoint[nZone][MAX_NUMBER_SLIDING][2];			/*!< \brief SlidingPoint[Periodic bc] and return the point that
                                                                         must be sent [0], and the matching point in the sliding bc[1]. */
	vector<unsigned long> SlidingElem[nZone][MAX_NUMBER_SLIDING];				/*!< \brief SlidingElem[Periodic bc] and return the elements that
                                                                       must be sent. */
  vector<unsigned long> NewBoundaryPoints[nZone][MAX_NUMBER_MARKER];  /*!< \brief Vector containing new points appearing on multiple boundaries. */
	vector<unsigned long> OldBoundaryElems[nZone][MAX_NUMBER_MARKER];  /*!< \brief Vector of old boundary elements. */

	long *PeriodicDomainIndex = NULL;
  bool **SlidingBC = NULL;
  
	/*--- Boolean vector to flag the points that belong to a sliding marker
   within each zone. This will be used in the logic for which points need
   to be replicated and added for the interfaces. ---*/
	SlidingBC = new bool*[nZone];
  for (iZone = 0; iZone < nZone; iZone++) {
    SlidingBC[iZone] = new bool[geometry_container[iZone]->GetnPoint()];
    
    /*--- Initialize all points to false ---*/
    for (iPoint = 0; iPoint < geometry_container[iZone]->GetnPoint(); iPoint++)
      SlidingBC[iZone][iPoint] = false;
    
    /*--- Set bool to true if point lies on a sliding interface ---*/
    for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++)
      if (config_container[iZone]->GetMarker_All_Boundary(iMarker) == SLIDING_INTERFACE)
          for (iVertex = 0; iVertex < geometry_container[iZone]->nVertex[iMarker]; iVertex++) {
            iPoint = geometry_container[iZone]->vertex[iMarker][iVertex]->GetNode();
            SlidingBC[iZone][iPoint] = true;
          }
	}
  
  /*--- Now that all of the sliding points have been tagged, determine the
   new points that will be added to each sliding boundary as a rind layer. ---*/
  for (iZone = 0; iZone < nZone; iZone++)
    for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++)
      if (config_container[iZone]->GetMarker_All_Boundary(iMarker) == SLIDING_INTERFACE) {
        
        /*--- Get internal index & donor zone for this sliding boundary ---*/
        iSliding = config_container[iZone]->GetMarker_All_Sliding(iMarker);
        jZone    = config_container[iZone]->GetSlideDonor_Zone(config_container[iZone]->GetMarker_All_Tag(iMarker));
        
        for (iVertex = 0; iVertex <  geometry_container[iZone]->nVertex[iMarker]; iVertex++) {
          
          /*--- iPoint is the original point on the surface and jPoint is the
           matching point on the sliding boundary in the adjacent zone.---*/
          iPoint = geometry_container[iZone]->vertex[iMarker][iVertex]->GetNode();
          jPoint = geometry_container[iZone]->vertex[iMarker][iVertex]->GetDonorPoint();
          
          /*--- Find & store any neighbor points to the sliding boundary in the donor zone (jZone). ---*/
          for (jNeighbor = 0; jNeighbor < geometry_container[jZone]->node[jPoint]->GetnPoint(); jNeighbor++) {
            kPoint = geometry_container[jZone]->node[jPoint]->GetPoint(jNeighbor);
            
            /*--- Only add points that don't already belong to the sliding boundary condition ---*/
            if (!SlidingBC[jZone][kPoint]) SlidingPoint[iZone][iSliding][0].push_back(kPoint);
          }
          
          /*--- Also store any elements that share a point with the sliding boundary ---*/
          for (jNeighbor = 0; jNeighbor < geometry_container[jZone]->node[jPoint]->GetnElem(); jNeighbor++) {
            jElem = geometry_container[jZone]->node[jPoint]->GetElem(jNeighbor);
            SlidingElem[iZone][iSliding].push_back(jElem);
          }
        }
      }


  /*--- Sort the points & elems and remove any duplicates by zone ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    
    /*--- Sort the points that must be sent by index and delete any repeats ---*/
    for (jMarker = 0; jMarker < config_container[iZone]->GetnMarker_All(); jMarker++)
      if (config_container[iZone]->GetMarker_All_Boundary(jMarker) == SLIDING_INTERFACE) {
        iSliding = config_container[iZone]->GetMarker_All_Sliding(jMarker);
        sort( SlidingPoint[iZone][iSliding][0].begin(), SlidingPoint[iZone][iSliding][0].end());
        IterPoint[iSliding][0] = unique( SlidingPoint[iZone][iSliding][0].begin(), SlidingPoint[iZone][iSliding][0].end());
        SlidingPoint[iZone][iSliding][0].resize( IterPoint[iSliding][0] - SlidingPoint[iZone][iSliding][0].begin() );
      }
    
    /*--- Create a list of the points that will receive information (new points).
     Note that the indexing for these new points starts immediately after 
     the last point in the original zone. ---*/
    iPointSliding = geometry_container[iZone]->GetnPoint();
    for (jMarker = 0; jMarker < config_container[iZone]->GetnMarker_All(); jMarker++)
      if (config_container[iZone]->GetMarker_All_Boundary(jMarker) == SLIDING_INTERFACE) {
        iSliding = config_container[iZone]->GetMarker_All_Sliding(jMarker);
        for (iPoint = 0; iPoint < SlidingPoint[iZone][iSliding][0].size(); iPoint++) {
          SlidingPoint[iZone][iSliding][1].push_back(iPointSliding);
          iPointSliding++;
        }
      }
    
    /*--- Sort and remove duplicates among the new elements too. ---*/
    for (jMarker = 0; jMarker < config_container[iZone]->GetnMarker_All(); jMarker++)
      if (config_container[iZone]->GetMarker_All_Boundary(jMarker) == SLIDING_INTERFACE) {
        iSliding = config_container[iZone]->GetMarker_All_Sliding(jMarker);
        sort( SlidingElem[iZone][iSliding].begin(), SlidingElem[iZone][iSliding].end());
        IterElem[iSliding] = unique( SlidingElem[iZone][iSliding].begin(), SlidingElem[iZone][iSliding].end());
        SlidingElem[iZone][iSliding].resize( IterElem[iSliding] - SlidingElem[iZone][iSliding].begin() );
      }
    
  }
  
  
//  /*--- Find any additional boundary elements that were left out of the
//   initial search because their points lie on multiple boundaries (corners
//   of the zone). ---*/
//	for (iZone = 0; iZone < nZone; iZone++) {
//    
//    /*--- Check all new points to see if they also lie on another boundary. ---*/
//    for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++)
//      for (iVertex = 0; iVertex < geometry_container[iZone]->nVertex[iMarker]; iVertex++) {
//        
//        /*--- iPoint is a node that lies on the current marker. ---*/
//        iPoint = geometry_container[iZone]->vertex[iMarker][iVertex]->GetNode();
//        
//        /*--- Check if iPoint lies on any other boundaries in this zone. ---*/
//        for (jMarker = 0; jMarker < config_container[iZone]->GetnMarker_All(); jMarker++)
//          if (config_container[iZone]->GetMarker_All_Boundary(jMarker) == SLIDING_INTERFACE) {
//            
//            /*--- An internal integer identifying this sliding boundary ---*/
//            iSliding = config_container[iZone]->GetMarker_All_Sliding(jMarker);
//            
//            /*--- jPoint is the SEND point. ---*/
//            for (kPoint = 0; kPoint < SlidingPoint[iZone][iSliding][0].size(); kPoint++) {
//              
//              jPoint = SlidingPoint[iZone][iSliding][0][kPoint];
//              
//              /*--- If the two match, then jPoint lies on this boundary.
//               However, we are concerned with the new points, so we
//               will store pPoint instead. ---*/
//              if (iPoint == jPoint) {
//                pPoint = SlidingPoint[iZone][iSliding][1][kPoint];
//                
//                /*--- We also want the type of boundary element that this point
//                 was within, so that we know what type of element to add
//                 built from the new points. ---*/
//                bool isJPoint, isSliding;
//                for(jElem = 0; jElem < geometry_container[iZone]->nElem_Bound[iMarker]; jElem++) {
//                  isJPoint = false; isSliding = false;
//                  
//                  //! Does this need to be in jZone?!?!
//                  for (iNode = 0; iNode < geometry_container[iZone]->bound[iMarker][jElem]->GetnNodes(); iNode++) {
//                    if (geometry_container[iZone]->bound[iMarker][jElem]->GetNode(iNode) == jPoint) isJPoint = true;
//                    if (SlidingBC[iZone][geometry_container[iZone]->bound[iMarker][jElem]->GetNode(iNode)]) isSliding = true;
//                  }
//                  
//                  /*--- If both points were found, store this element. ---*/
//                  if (isJPoint && isSliding) OldBoundaryElems[iZone][iMarker].push_back(jElem);
//                  
//                }
//              }
//              
//            }
//          }
//      }
//  }
  
  /*--- Sort & remove any duplicates among this new set of additional elements. ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++) {
      sort( OldBoundaryElems[iZone][iMarker].begin(), OldBoundaryElems[iZone][iMarker].end());
      IterNewElem[iMarker] = unique( OldBoundaryElems[iZone][iMarker].begin(), OldBoundaryElems[iZone][iMarker].end());
      OldBoundaryElems[iZone][iMarker].resize( IterNewElem[iMarker] - OldBoundaryElems[iZone][iMarker].begin());
    }
  }
  
  
  /*--- Create the new boundary elements. Points making up these new
   elements must either be SEND/RECEIVE or sliding points. ---*/
  
  CPrimalGrid**** newBound = new CPrimalGrid***[nZone];            /*!< \brief Boundary vector for new periodic elements (primal grid information). */
	unsigned long **nNewElem_Bound = new unsigned long*[nZone];			/*!< \brief Number of new periodic elements of the boundary. */
  
  for (iZone = 0; iZone < nZone; iZone++) {
    
    nNewElem_Bound[iZone] = new unsigned long[config_container[iZone]->GetnMarker_All()];
    newBound[iZone] = new CPrimalGrid**[config_container[iZone]->GetnMarker_All()];
    
    for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++) {
        
        nNewElem_Bound[iZone][iMarker] = OldBoundaryElems[iZone][iMarker].size();
        newBound[iZone][iMarker]       = new CPrimalGrid*[nNewElem_Bound[iZone][iMarker]];
        
        /*--- Loop through all new elements to be added. ---*/
        for (iElem = 0; iElem < nNewElem_Bound[iZone][iMarker]; iElem++) {
          jElem = OldBoundaryElems[iZone][iMarker][iElem];
          
          /*--- Loop through all nodes of this element. ---*/
          for (iNode = 0; iNode < geometry_container[iZone]->bound[iMarker][jElem]->GetnNodes(); iNode++) {
            pPoint = geometry_container[iZone]->bound[iMarker][jElem]->GetNode(iNode);
            
            /*--- Check if this node is a send point. If so, the corresponding
             receive point will be used in making the new boundary element. ---*/
            for (jMarker = 0; jMarker < config_container[iZone]->GetnMarker_All(); jMarker++)
              if (config_container[iZone]->GetMarker_All_Boundary(jMarker) == SLIDING_INTERFACE) {
                  
                  /*--- An integer identify the sliding boundary condition ---*/
                  iSliding = config_container[iZone]->GetMarker_All_Sliding(jMarker);
                  for (kElem = 0; kElem < SlidingPoint[iZone][iSliding][0].size(); kElem++) {
                    if (pPoint == SlidingPoint[iZone][iSliding][0][kElem])
                      newNodes[iNode] = SlidingPoint[iZone][iSliding][1][kElem];
                  }
                }
            
            /*--- Check if this node is a sliding point. If so, the corresponding
             sliding point will be used in making the new boundary element. ---*/
            if (SlidingBC[iZone][pPoint]) {
              
              /*--- Find the corresponding sliding point. ---*/
              for (jMarker = 0; jMarker < config_container[iZone]->GetnMarker_All(); jMarker++) {
                if (config_container[iZone]->GetMarker_All_Boundary(jMarker) == SLIDING_INTERFACE) {
                    
                    for (iVertex = 0; iVertex < geometry_container[iZone]->nVertex[jMarker]; iVertex++) {
                      if (pPoint == geometry_container[iZone]->vertex[jMarker][iVertex]->GetNode()) {
                        kMarker = jMarker;
                        jVertex = iVertex;
                      }
                    }
                  }
              }
              newNodes[iNode] = geometry_container[iZone]->vertex[kMarker][jVertex]->GetDonorPoint();
            }
          }
          
          /*--- Now instantiate the new element. ---*/
          VTK_Type = geometry_container[iZone]->bound[iMarker][jElem]->GetVTK_Type();
          switch(VTK_Type) {
            case LINE:
              newBound[iZone][iMarker][iElem] = new CLine(newNodes[0],newNodes[1],2);
              break;
            case TRIANGLE:
              newBound[iZone][iMarker][iElem] = new CTriangle(newNodes[0],newNodes[1],newNodes[2],3);
              break;
            case RECTANGLE:
              newBound[iZone][iMarker][iElem] = new CRectangle(newNodes[0],newNodes[1],newNodes[2],newNodes[3],3);
              break;
          }
          
        }
      }
    
    
  }
  
  //	delete [] SlidingBC;
  
  //	/*--- Set periodic boundary conditions ---*/
  //	geometry->SetPeriodicBoundary(config);
  //
//  /*--- Original grid in .vtk format for debugigng purposes ---*/
//  for (iZone = 0; iZone < nZone; iZone++) {
//
//    //  strcpy (cstr, config_container[iZone]->GetSubDomain_FileName(iZone).c_str());
//    char buffer[100];
//    sprintf (buffer, "domain_orig_%d.plt", int(iZone));
//    //  strcat(cstr,buffer);
//    //  strcpy (buffer_vtk, "periodic_orig.vtk");
//    geometry_container[iZone]->SetTecPlot(buffer);
//  }
  
  
//#ifdef debug
  //////////
  /*---- BUILDING NEW MESH FILES STARTS HERE ---*/
  //////////
  unsigned short nSliding = 0;
  unsigned long nelem_pyramid, nelem_quad, nelem_tetra, nelem_hexa, nelem_triangle, nelem_wedge;
  unsigned long newElementsBound, nElem_new, nPoint_new, nElem_Storage;
  CPrimalGrid*** elem;	/*!< \brief Element vector (primal grid information). */
	CPrimalGrid**** bound;	/*!< \brief Boundary vector (primal grid information). */
	CPoint*** node;			/*!< \brief Node vector (dual grid information). */
	CEdge** edge;			/*!< \brief Edge vector (dual grid information). */
	CVertex*** vertex;		/*!< \brief Boundary Vertex vector (dual grid information). */
	unsigned long *nVertex;	/*!< \brief Number of vertex for each marker. */
  unsigned long **Index;
  unsigned long **nElem_Bound;
  unsigned long *nElem = new unsigned long[nZone];
  unsigned long *nPoint = new unsigned long[nZone];
  elem = new CPrimalGrid**[nZone];
  bound = new CPrimalGrid***[nZone];
  node = new CPoint**[nZone];
  Index = new unsigned long*[nZone];
  nElem_Bound = new unsigned long*[nZone];
  string *Tag_to_Marker = new string [MAX_INDEX_VALUE];
  double *Coord_j;
  
  /*--- Find max number of points in any zone ---*/
  unsigned long Max_Points = 0;
  for (iZone = 0; iZone < nZone; iZone++) {
    if (geometry_container[iZone]->GetnPoint() > Max_Points)
      Max_Points = geometry_container[iZone]->GetnPoint();
  }
  
  
  for (iZone = 0; iZone < nZone; iZone++) {
    
    /*--- Compute the number of sliding boundaries in this zone ---*/
    for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++)
      if (config_container[iZone]->GetMarker_All_Boundary(iMarker) == SLIDING_INTERFACE)
        nSliding++;
    
    /*--- Write the number of dimensions of the problem ---*/
    nDim = geometry_container[iZone]->GetnDim();
    
    //  /*--- Copy the new boundary element information from the geometry class.
    //   Be careful, as these are pointers to vectors/objects. ---*/
    //  nNewElem_BoundPer = geometry->nNewElem_Bound;
    //  newBoundPer       = geometry->newBound;
    
    /*--- Count the total number of new boundary elements in this zone. ---*/
    newElementsBound = 0;
    for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++)
      newElementsBound += nNewElem_Bound[iZone][iMarker];
    
    /*--- Loop over the original grid to perform the dimensionalizaton of the new vectors ---*/
    nElem_new = 0; nPoint_new = 0;
    for (jMarker = 0; jMarker < config_container[iZone]->GetnMarker_All(); jMarker++)
      if (config_container[iZone]->GetMarker_All_Boundary(jMarker) == SLIDING_INTERFACE) {
        
        /*--- An integer identifying the sliding boundary condition ---*/
        iSliding = config_container[iZone]->GetMarker_All_Sliding(jMarker);
          nElem_new  += SlidingElem[iZone][iSliding].size();
          nPoint_new += SlidingPoint[iZone][iSliding][0].size();
        
        /*--- Print some info about the boundaries to the console. ---*/
        cout << "Constructing new geometry for Zone " << iZone+1 << "." << endl;
        cout << "Number of new points: " << nPoint_new << "." << endl;
        cout << "Number of new interior elements: " << nElem_new << "." << endl;
        cout << "Number of new boundary elements added to preexisting markers: " << newElementsBound << "." << endl;
        
      }
    
    /*--- Create a copy of the original grid ---*/
    elem[iZone] = new CPrimalGrid*[geometry_container[iZone]->GetnElem() + nElem_new];
    nelem_pyramid = 0; nelem_quad = 0; nelem_tetra = 0; nelem_hexa = 0; nelem_triangle = 0; nelem_wedge = 0;
    
    for (iElem = 0; iElem < geometry_container[iZone]->GetnElem(); iElem ++) {
      switch(geometry_container[iZone]->elem[iElem]->GetVTK_Type()) {
        case TRIANGLE:
          elem[iZone][iElem] = new CTriangle(geometry_container[iZone]->elem[iElem]->GetNode(0),
                                             geometry_container[iZone]->elem[iElem]->GetNode(1),
                                             geometry_container[iZone]->elem[iElem]->GetNode(2), 2);
          nelem_triangle++;
          break;
          
        case RECTANGLE:
          elem[iZone][iElem] = new CRectangle(geometry_container[iZone]->elem[iElem]->GetNode(0),
                                              geometry_container[iZone]->elem[iElem]->GetNode(1),
                                              geometry_container[iZone]->elem[iElem]->GetNode(2),
                                              geometry_container[iZone]->elem[iElem]->GetNode(3), 2);
          nelem_quad++;
          break;
          
        case TETRAHEDRON:
          elem[iZone][iElem] = new CTetrahedron(geometry_container[iZone]->elem[iElem]->GetNode(0),
                                                geometry_container[iZone]->elem[iElem]->GetNode(1),
                                                geometry_container[iZone]->elem[iElem]->GetNode(2),
                                                geometry_container[iZone]->elem[iElem]->GetNode(3));
          nelem_tetra++;
          break;
          
        case HEXAHEDRON:
          elem[iZone][iElem] = new CHexahedron(geometry_container[iZone]->elem[iElem]->GetNode(0),
                                               geometry_container[iZone]->elem[iElem]->GetNode(1),
                                               geometry_container[iZone]->elem[iElem]->GetNode(2),
                                               geometry_container[iZone]->elem[iElem]->GetNode(3),
                                               geometry_container[iZone]->elem[iElem]->GetNode(4),
                                               geometry_container[iZone]->elem[iElem]->GetNode(5),
                                               geometry_container[iZone]->elem[iElem]->GetNode(6),
                                               geometry_container[iZone]->elem[iElem]->GetNode(7));
          nelem_hexa++;
          break;
          
        case WEDGE:
          elem[iZone][iElem] = new CWedge(geometry_container[iZone]->elem[iElem]->GetNode(0),
                                          geometry_container[iZone]->elem[iElem]->GetNode(1),
                                          geometry_container[iZone]->elem[iElem]->GetNode(2),
                                          geometry_container[iZone]->elem[iElem]->GetNode(3),
                                          geometry_container[iZone]->elem[iElem]->GetNode(4),
                                          geometry_container[iZone]->elem[iElem]->GetNode(5));
          nelem_wedge++;
          break;
          
        case PYRAMID:
          elem[iZone][iElem] = new CPyramid(geometry_container[iZone]->elem[iElem]->GetNode(0),
                                            geometry_container[iZone]->elem[iElem]->GetNode(1),
                                            geometry_container[iZone]->elem[iElem]->GetNode(2),
                                            geometry_container[iZone]->elem[iElem]->GetNode(3),
                                            geometry_container[iZone]->elem[iElem]->GetNode(4));
          nelem_pyramid++;
          break;
          
      }
    }
    
    /*--- Create a list with all the points and the new index values ---*/
    // needs to have enough space for ALL zones
    
    Index[iZone] = new unsigned long[Max_Points];
    for (iPoint = 0; iPoint < Max_Points; iPoint ++)
      Index[iZone][iPoint] = 0;
    
    for (jMarker = 0; jMarker < config_container[iZone]->GetnMarker_All(); jMarker++)
      if (config_container[iZone]->GetMarker_All_Boundary(jMarker) == SLIDING_INTERFACE) {
        
        /*--- An integer identifying the sliding boundary condition ---*/
        iSliding = config_container[iZone]->GetMarker_All_Sliding(jMarker);
        
        for (iNeighbor = 0; iNeighbor < SlidingPoint[iZone][iSliding][0].size(); iNeighbor++) {
          //! same flip as below
          iPoint =  SlidingPoint[iZone][iSliding][0][iNeighbor];
          Index[iZone][iPoint] = SlidingPoint[iZone][iSliding][1][iNeighbor];
        }
      }
    
    for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++)
      if (config_container[iZone]->GetMarker_All_Boundary(iMarker) == SLIDING_INTERFACE) {
        for (iVertex = 0; iVertex < geometry_container[iZone]->GetnVertex(iMarker); iVertex++) {
          //! Note that iPoint/jPoint are flipped from the periodic condition
          //! because now we aren't working within the original domain. We have
          //! to include the original nodes on the sliding surface in iZone.
          jPoint = geometry_container[iZone]->vertex[iMarker][iVertex]->GetNode();
          iPoint = geometry_container[iZone]->vertex[iMarker][iVertex]->GetDonorPoint();
          Index[iZone][iPoint] = jPoint;
        }
      }
    
    /*--- Add the new elements due to the sliding boundary condtion ---*/
    iElem = geometry_container[iZone]->GetnElem();
    
    for (jMarker = 0; jMarker < config_container[iZone]->GetnMarker_All(); jMarker++)
      if (config_container[iZone]->GetMarker_All_Boundary(jMarker) == SLIDING_INTERFACE) {
        
        /*--- An integer identify the periodic boundary condition ---*/
        iSliding = config_container[iZone]->GetMarker_All_Sliding(jMarker);
        
        /*--- Also need to pull information from the correct domain ---*/
        jZone = config_container[iZone]->GetSlideDonor_Zone(config_container[iZone]->GetMarker_All_Tag(jMarker));
        
        for (iNeighbor = 0; iNeighbor < SlidingElem[iZone][iSliding].size(); iNeighbor++) {
          jElem = SlidingElem[iZone][iSliding][iNeighbor];
          
          switch(geometry_container[jZone]->elem[jElem]->GetVTK_Type()) {
            case TRIANGLE:
              elem[iZone][iElem] = new CTriangle(Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(0)],
                                                 Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(1)],
                                                 Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(2)], 2);
              iElem++; nelem_triangle++;
              break;
              
            case RECTANGLE:
              elem[iZone][iElem] = new CRectangle(Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(0)],
                                                  Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(1)],
                                                  Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(2)],
                                                  Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(3)], 2);
              iElem++; nelem_quad++;
              break;
              
            case TETRAHEDRON:
              elem[iZone][iElem] = new CTetrahedron(Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(0)],
                                                    Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(1)],
                                                    Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(2)],
                                                    Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(3)]);
              iElem++; nelem_tetra++;
              break;
              
            case HEXAHEDRON:
              elem[iZone][iElem] = new CHexahedron(Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(0)],
                                                   Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(1)],
                                                   Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(2)],
                                                   Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(3)],
                                                   Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(4)],
                                                   Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(5)],
                                                   Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(6)],
                                                   Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(7)]);
              iElem++; nelem_hexa++;
              break;
              
            case WEDGE:
              elem[iZone][iElem] = new CWedge(Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(0)],
                                              Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(1)],
                                              Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(2)],
                                              Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(3)],
                                              Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(4)],
                                              Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(5)]);
              iElem++; nelem_wedge++;
              break;
              
            case PYRAMID:
              elem[iZone][iElem] = new CPyramid(Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(0)],
                                                Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(1)],
                                                Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(2)],
                                                Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(3)],
                                                Index[iZone][geometry_container[jZone]->elem[jElem]->GetNode(4)]);
              iElem++; nelem_pyramid++;
              break;
              
          }
        }
      }
    
    
    nElem_Storage = nelem_triangle*4 + nelem_quad*5 + nelem_tetra*5 + nelem_hexa*9 + nelem_wedge*7 + nelem_pyramid*6;
    nElem[iZone] = geometry_container[iZone]->GetnElem() + nElem_new;
    
    /*--- Add the old points ---*/
    node[iZone] = new CPoint*[geometry_container[iZone]->GetnPoint() + nPoint_new];
    for (iPoint = 0; iPoint < geometry_container[iZone]->GetnPoint(); iPoint++) {
      if (geometry_container[iZone]->GetnDim() == 2)
        node[iZone][iPoint] = new CPoint(geometry_container[iZone]->node[iPoint]->GetCoord(0),
                                         geometry_container[iZone]->node[iPoint]->GetCoord(1), iPoint, config);
      if (geometry_container[iZone]->GetnDim() == 3)
        node[iZone][iPoint] = new CPoint(geometry_container[iZone]->node[iPoint]->GetCoord(0),
                                         geometry_container[iZone]->node[iPoint]->GetCoord(1),
                                         geometry_container[iZone]->node[iPoint]->GetCoord(2), iPoint, config);
    }
    
    /*--- Add the new points due to the sliding boundary condtion ---*/
    for (jMarker = 0; jMarker < config_container[iZone]->GetnMarker_All(); jMarker++)
      if (config_container[iZone]->GetMarker_All_Boundary(jMarker) == SLIDING_INTERFACE) {
        
        /*--- An integer identify the periodic boundary condition ---*/
        iSliding = config_container[iZone]->GetMarker_All_Sliding(jMarker);
        
        /*--- Also need to pull information from the correct domain ---*/
        jZone = config_container[iZone]->GetSlideDonor_Zone(config_container[iZone]->GetMarker_All_Tag(jMarker));
        
        for (iNeighbor = 0; iNeighbor < SlidingPoint[iZone][iSliding][0].size(); iNeighbor++) {
          
          //      cout << iPoint << endl;
          /*--- Retrieve node information for this boundary point. ---*/
          jPoint = SlidingPoint[iZone][iSliding][0][iNeighbor];
          iPoint = SlidingPoint[iZone][iSliding][1][iNeighbor];
          //			Coord_i = geometry[iZone]->node[iPoint]->GetCoord();
          //cout << iPoint << " " <<jPoint << endl;
          
          Coord_j = geometry_container[jZone]->node[jPoint]->GetCoord();
          
          /////////
          ///////// Double check this later, make sure i'm using the correct domains/points
          /////////
          //cout << Coord_j[0] << endl;
          //cout << iPoint << " " <<jPoint << endl;
          
          /*--- Save the new points with the new coordinates. ---*/
          if (geometry_container[iZone]->GetnDim() == 2)
            node[iZone][iPoint] = new CPoint(Coord_j[0], Coord_j[1], iPoint, config);
          if (geometry_container[iZone]->GetnDim() == 3)
            node[iZone][iPoint] = new CPoint(Coord_j[0], Coord_j[1], Coord_j[2], iPoint, config);
          
        }
      }
    nPoint[iZone] = geometry_container[iZone]->GetnPoint() + nPoint_new;
    
    
    /*--- Add the old boundary, reserving space for two new bc (send/recive sliding bc) ---*/
    unsigned short nMarker = config_container[iZone]->GetnMarker_All();
//    for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++)
//        nMarker++;
//      }
    nMarker += 2;
    nElem_Bound[iZone] = new unsigned long [nMarker];
    bound[iZone] = new CPrimalGrid**[nMarker];
    
    /////// Setting new # of markers in config
    //	config_container[iZone]->SetnMarker_All(nMarker);
    
    /*--- Copy the old boundary --- note we're restarting indexing from  0 ---*/
    jMarker = 0;
    for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++) {
      
      bound[iZone][jMarker] = new CPrimalGrid* [geometry_container[iZone]->GetnElem_Bound(iMarker)];
      
      for (iVertex = 0; iVertex < geometry_container[iZone]->GetnElem_Bound(iMarker); iVertex++) {
        
        if (geometry_container[iZone]->bound[iMarker][iVertex]->GetVTK_Type() == LINE)
          bound[iZone][jMarker][iVertex] = new CLine(geometry_container[iZone]->bound[iMarker][iVertex]->GetNode(0),
                                                     geometry_container[iZone]->bound[iMarker][iVertex]->GetNode(1), 2);
        
        if (geometry_container[iZone]->bound[iMarker][iVertex]->GetVTK_Type() == TRIANGLE)
          bound[iZone][jMarker][iVertex] = new CTriangle(geometry_container[iZone]->bound[iMarker][iVertex]->GetNode(0),
                                                         geometry_container[iZone]->bound[iMarker][iVertex]->GetNode(1),
                                                         geometry_container[iZone]->bound[iMarker][iVertex]->GetNode(2), 3);
        
        if (geometry_container[iZone]->bound[iMarker][iVertex]->GetVTK_Type() == RECTANGLE)
          bound[iZone][jMarker][iVertex] = new CRectangle(geometry_container[iZone]->bound[iMarker][iVertex]->GetNode(0),
                                                          geometry_container[iZone]->bound[iMarker][iVertex]->GetNode(1),
                                                          geometry_container[iZone]->bound[iMarker][iVertex]->GetNode(2),
                                                          geometry_container[iZone]->bound[iMarker][iVertex]->GetNode(3), 3);
        
      }
      
      nElem_Bound[iZone][jMarker] = geometry_container[iZone]->GetnElem_Bound(iMarker);
      Tag_to_Marker[jMarker] = geometry_container[iZone]->GetMarker_Tag(iMarker);
      jMarker++;
    }
  
cout << jMarker << endl;
cout << nMarker << endl;
  
    /////////
    ///////// Setting SEND/RECEIVE - this will change in the future to something else
    ////////  to label the points for interpolation
    
    unsigned short iMarkerSend, iMarkerReceive, jSliding;
    unsigned long  Counter_Send = 0, Counter_Receive = 0;
    string iTag, jTag;
    
    /*--- First compute the Send/Receive boundaries, count the number of points ---*/
    Counter_Send = 0; 	Counter_Receive = 0;
    for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++)
      if (config_container[iZone]->GetMarker_All_Boundary(iMarker) == SLIDING_INTERFACE) {
        
        cout << "counting points for SEND/RECV" <<endl;
        
        /*--- Tag & internal index of the current sliding boundary marker ---*/
        iTag     = config_container[iZone]->GetMarker_All_Tag(iMarker);
        iSliding = config_container[iZone]->GetMarker_All_Sliding(iMarker);
        
        /*--- Index values for RECEIVE points are all within iZone ---*/
        if (SlidingPoint[iZone][iSliding][1].size() != 0)
          Counter_Receive += SlidingPoint[iZone][iSliding][1].size();
        
        /*--- Get zone number, tag, & internal index for the donor marker. ---*/
        jZone    = config_container[iZone]->GetSlideDonor_Zone(iTag);
        jTag     = config_container[iZone]->GetMarker_Sliding_Donor(iTag);
        jMarker  = config_container[jZone]->GetTag_Marker_All(jTag);
        jSliding = config_container[jZone]->GetMarker_All_Sliding(jMarker);

        /*--- Index values from SEND points in iZone are stored under jZone/jSliding. ---*/
        if (SlidingPoint[jZone][jSliding][0].size() != 0)
          Counter_Send += SlidingPoint[jZone][jSliding][0].size();

        
        cout << "Send: "<<Counter_Send <<" RECV: "<< Counter_Receive <<endl;
      }
    
    /*--- Adimensionalization of the new boundaries ---*/
    iMarkerSend = nMarker - 2; iMarkerReceive = nMarker - 1;
    config_container[iZone]->SetMarker_All_SendRecv(iMarkerSend,1);
    config_container[iZone]->SetMarker_All_SendRecv(iMarkerReceive,-1);
    nElem_Bound[iZone][iMarkerSend] = Counter_Send;
    nElem_Bound[iZone][iMarkerReceive] = Counter_Receive;
    bound[iZone][iMarkerSend] = new CPrimalGrid* [Counter_Send];
    bound[iZone][iMarkerReceive] = new CPrimalGrid* [Counter_Receive];
    
    cout << "set new bound sizes" <<endl;
    /*--- First we do the send ---*/
    iVertex = 0;
    for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++)
      if (config_container[iZone]->GetMarker_All_Boundary(iMarker) == SLIDING_INTERFACE) {
        
        /*--- Tag & internal index of the current sliding boundary marker ---*/
        iTag     = config_container[iZone]->GetMarker_All_Tag(iMarker);
        iSliding = config_container[iZone]->GetMarker_All_Sliding(iMarker);
        
        /*--- Get zone number, tag, & internal index for the donor marker. ---*/
        jZone    = config_container[iZone]->GetSlideDonor_Zone(iTag);
        jTag     = config_container[iZone]->GetMarker_Sliding_Donor(iTag);
        jMarker  = config_container[jZone]->GetTag_Marker_All(jTag);
        jSliding = config_container[jZone]->GetMarker_All_Sliding(jMarker);
                
        if (SlidingPoint[jZone][jSliding][0].size() != 0)
          for (iNeighbor = 0; iNeighbor < SlidingPoint[jZone][jSliding][0].size(); iNeighbor++) {
            bound[iZone][iMarkerSend][iVertex] = new CVertexMPI(SlidingPoint[jZone][jSliding][0][iNeighbor], 3);
            bound[iZone][iMarkerSend][iVertex]->SetRotation_Type(0);
            bound[iZone][iMarkerSend][iVertex]->SetMatching_Zone(jZone);
            iVertex++;
          }
      }
    
    /*--- Second we do the receive ---*/
    iVertex = 0;
    for (jMarker = 0; jMarker < config_container[iZone]->GetnMarker_All(); jMarker++)
      if (config_container[iZone]->GetMarker_All_Boundary(jMarker) == SLIDING_INTERFACE) {
        
        /*--- An integer identifying the periodic boundary condition ---*/
        iSliding = config_container[iZone]->GetMarker_All_Sliding(jMarker);
        jZone = config_container[iZone]->GetSlideDonor_Zone(config_container[iZone]->GetMarker_All_Tag(jMarker));
        
        if (SlidingPoint[iZone][iSliding][1].size() != 0)
          for (iNeighbor = 0; iNeighbor < SlidingPoint[iZone][iSliding][1].size(); iNeighbor++) {
            bound[iZone][iMarkerReceive][iVertex] = new CVertexMPI(SlidingPoint[iZone][iSliding][1][iNeighbor], 3);
            bound[iZone][iMarkerReceive][iVertex]->SetRotation_Type(0);
            bound[iZone][iMarkerReceive][iVertex]->SetMatching_Zone(jZone);
            iVertex++;
          }
      }
    
  }
    
  /*--- Write tecplot files for visualizing the sliding interfaces ---*/
  SetTecplot(node, nPoint, elem, nElem, config_container, nZone, nDim);
  
  /*--- Write the SU2 mesh file with the new sliding interfaces ---*/
  SetMeshFile(node, nPoint, elem, nElem, bound, nElem_Bound, newBound,
              nNewElem_Bound, config_container, nZone, nDim);
	
	return 1;
}


unsigned short GetnZone(string val_mesh_filename, unsigned short val_format, CConfig *config) {
  
	string text_line, Marker_Tag;
	ifstream mesh_file;
	short nZone = 1;
	bool isFound = false;
	char cstr[200];
	string::size_type position;
	int rank = MASTER_NODE;

  
	switch (val_format) {
    case SU2:
      
      /*--- Open grid file ---*/
      strcpy (cstr, val_mesh_filename.c_str());
      mesh_file.open(cstr, ios::in);
      if (mesh_file.fail()) {
        cout << "There is no geometry file (GetnZone))!" << endl;
        cout << "Press any key to exit..." << endl;
        cin.get();
        exit(1);
      }
      
      /*--- Open the SU2 mesh file ---*/
      while (getline (mesh_file,text_line)) {
        
        /*--- Search for the "NZONE" keyword to see if there are multiple Zones ---*/
        position = text_line.find ("NZONE=",0);
        if (position != string::npos) {
          text_line.erase (0,6); nZone = atoi(text_line.c_str()); isFound = true;
          if (rank == MASTER_NODE) {
            if (nZone <= 1) {
              cout << "Error: Number of mesh zones is less than 1 !!!" << endl;
              cout << "Press any key to exit..." << endl;
              cin.get();
              exit(1);
            }
          }
        }
      }
      
      /*--- If the "NZONE" keyword was not found, throw an error. We need
       more than one zone in order to build the sliding interface. ---*/
      if (!isFound) {
        if (rank == MASTER_NODE) {
          cout << "Error: Number of mesh zones is not greater than 1 !!!" << endl;
          cout << "Press any key to exit..." << endl;
          cin.get();
          exit(1);
        }
      }
      break;
      
    case CGNS:
      cout << "Sliding mesh capability currently unavailable with CGNS!" << endl;
      cout << "Press any key to exit..." << endl;
      cin.get();
      exit(1);
      break;
      
    case NETCDF_ASCII:
      cout << "Sliding mesh capability currently unavailable with NETCDF_ASCII!" << endl;
      cout << "Press any key to exit..." << endl;
      cin.get();
      exit(1);
      break;
      
	}
  
	return (unsigned short) nZone;
}

unsigned short GetnDim(string val_mesh_filename, unsigned short val_format) {
  
	string text_line, Marker_Tag;
	ifstream mesh_file;
	short nDim = 3;
	bool isFound = false;
	char cstr[200];
	string::size_type position;
	
	switch (val_format) {
    case SU2:
      
      /*--- Open grid file ---*/
      strcpy (cstr, val_mesh_filename.c_str());
      mesh_file.open(cstr, ios::in);
      
      /*--- Read SU2 mesh file ---*/
      while (getline (mesh_file,text_line)) {
        /*--- Search for the "NDIM" keyword to see if there are multiple Zones ---*/
        position = text_line.find ("NDIME=",0);
        if (position != string::npos) {
          text_line.erase (0,6); nDim = atoi(text_line.c_str()); isFound = true;
        }
      }
      break;
      
    case CGNS:
      nDim = 3;
      break;
      
    case NETCDF_ASCII:
      nDim = 3;
      break;
	}
	return (unsigned short) nDim;
}

void Geometric_Definition(CGeometry *geometry_container, CConfig *config_container) {
  
  /*--- Compute elements surrounding points, points surrounding points,
   and elements surrounding elements ---*/
  cout << "Setting local point and element connectivity." <<endl;
  geometry_container->SetEsuP();
  geometry_container->SetPsuP();
  geometry_container->SetEsuE();
  
  /*--- Check the orientation before computing geometric quantities ---*/
  cout << "Checking the numerical grid orientation." <<endl;
  geometry_container->SetBoundVolume();
  geometry_container->Check_Orientation(config_container);
  
  /*--- Create the edge structure ---*/
  cout << "Identifying edges and vertices." <<endl;
  geometry_container->SetEdges();
  geometry_container->SetVertex(config_container);
  
  /*--- Compute center of gravity ---*/
  cout << "Computing centers of gravity." << endl;
  geometry_container->SetCG();
  
  /*--- Create the control volume structures ---*/
  cout << "Setting the control volume structure." << endl;
  geometry_container->SetControlVolume(config_container, ALLOCATE);
  geometry_container->SetBoundControlVolume(config_container, ALLOCATE);
  
}

void MatchSliding_Interfaces(CGeometry **geometry_container, CConfig **config_container, unsigned short nZone) {
  
  /*--- Local variables ---*/
  bool isBadMatch = false;
  
  unsigned short iZone, jZone, iMarker, jMarker, kMarker = 0;
  unsigned short iSliding, iDim;
  
  unsigned long iPoint, jPoint, kPoint, pPoint;
  unsigned long iElem, jElem, kElem, iVertex, jVertex;
  
  double epsilon = 1e-10, mindist, dist, *Coord_i, *Coord_j;
  
  string iTag, jTag;
  
  /*--- Send an initial message to the console. ---*/
  cout << "Preprocessing the sliding mesh interfaces:\n" <<endl;
  
  /*--- Loop over all zones. ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    
  /*--- Loop through each marker to find any sliding boundaries. ---*/
  for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++)
    if (config_container[iZone]->GetMarker_All_Boundary(iMarker) == SLIDING_INTERFACE) {

      /*--- Tag of the current sliding boundary marker ---*/
      iTag  = config_container[iZone]->GetMarker_All_Tag(iMarker);
      
      /*--- Get zone number and tag for the donor marker. ---*/
      jZone = config_container[iZone]->GetSlideDonor_Zone(iTag);
      jTag  = config_container[iZone]->GetMarker_Sliding_Donor(iTag);
      
      /*--- Get the donor marker's index within the other zone ---*/
      jMarker = config_container[jZone]->GetTag_Marker_All(jTag);
            
      /*--- Write some info to the console. ---*/
      cout << "Matching marker '" << iTag;
      cout << "' in Zone " << iZone << " with donor '";
      cout << jTag << "' in Zone " << jZone << "." << endl;
      
      /*--- Loop through all vertices on this marker. ---*/
      for (iVertex = 0; iVertex < geometry_container[iZone]->nVertex[iMarker]; iVertex++) {
        
        /*--- Retrieve node information for this boundary point. ---*/
        iPoint  = geometry_container[iZone]->vertex[iMarker][iVertex]->GetNode();
        Coord_i = geometry_container[iZone]->node[iPoint]->GetCoord();
        
        /*--- Perform a brute-force search along the donor boundary
         marker in order to find the closest matching donor point. ---*/
        mindist = 1e10;
        for (jVertex = 0; jVertex < geometry_container[jZone]->nVertex[jMarker]; jVertex++) {
          
          /*--- Retrieve information for this jPoint on donor marker. ---*/
          jPoint  = geometry_container[jZone]->vertex[jMarker][jVertex]->GetNode();
          Coord_j = geometry_container[jZone]->node[jPoint]->GetCoord();
          
          /*--- Compute distance between iPoint & jPoint. ---*/
          dist = 0.0;
          for (iDim = 0; iDim < geometry_container[iZone]->GetnDim(); iDim++){
            dist += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
          }
          dist = sqrt(dist);
          
          /*---  Store node index if this is the closest point thus far. ---*/
          if (dist < mindist) { mindist = dist; pPoint = jPoint; }
        }
        
        /*--- Set the matching point for this iPoint. ---*/
        geometry_container[iZone]->vertex[iMarker][iVertex]->SetDonorPoint(pPoint);
        
        /*--- Print warning if the nearest point was not within
         the specified tolerance. Computation will continue. ---*/
        if (mindist > epsilon) {
          isBadMatch = true;
          cout.precision(10);
          cout << endl;
          cout << "   Non-match for point " << iPoint << ".\tNearest";
          cout << " donor distance: " << scientific << mindist << ".";
        }
      }
      
      /*--- Print final message in case there were any bad matches. ---*/
      if (isBadMatch) {
        cout << endl;
        cout << "\n !!! Warning !!!" << endl;
        cout << "Initial sliding interface not 1-to-1. Computation will continue.\n";
      } else {
        cout << "Initial sliding interface is 1-to-1 within specififed tolerance.\n";
      }
      cout << endl;
      isBadMatch = false;
      
    }
  }
}

void SetTecplot(CPoint ***node, unsigned long *nPoint, CPrimalGrid ***elem,
                unsigned long *nElem, CConfig **config_container, unsigned short nZone, unsigned short nDim) {
  
  /*--- Local variables ---*/
  unsigned short iDim, iZone;
  unsigned long iPoint, iElem;
  ofstream tecplot_file;
  char buffer[100];
  string filename;
  filename.assign(config_container[ZONE_0]->GetMesh_Out_FileName());
  filename.erase(filename.end()-4, filename.end());
  filename.append(".plt");
  
  strcpy (buffer, filename.c_str());
    tecplot_file.open(buffer, ios::out);
  
  /*--- Set the header once for a multi-zone Tecplot file. ---*/
  if (nDim == 2) {
    tecplot_file << " TITLE = \"Multi-zone mesh created by SU2_SMC\" " << endl;
    tecplot_file << " VARIABLES = \"x\",\"y\" " << endl;
  }
  if (nDim == 3) {
    tecplot_file << " TITLE = \"Multi-zone mesh created by SU2_SMC\" " << endl;
    tecplot_file << " VARIABLES = \"x\",\"y\",\"z\" " << endl;
  }
  
  /*--- Loop over all zones. ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    
    /*--- Header for each individual zone. ---*/
    if (nDim == 2) {
      tecplot_file << " ZONE T = \"Time = 0.0\", N= "<< nPoint[iZone] <<" , E = ";
      tecplot_file << nElem[iZone] <<" , F = FEPOINT, ET = QUADRILATERAL"<< endl;
    }
    if (nDim == 3) {
      tecplot_file << " ZONE T = \"Time = 0.0\", N= "<< nPoint[iZone] <<" , E = ";
      tecplot_file << nElem[iZone] <<" , F = FEPOINT, ET = BRICK"<< endl;
    }
    
    /*--- Write all point coordinates (x,y,z) first. ---*/
    for(iPoint = 0; iPoint < nPoint[iZone]; iPoint++) {
      for(iDim = 0; iDim < nDim; iDim++)
        tecplot_file << scientific << node[iZone][iPoint]->GetCoord(iDim) << "\t";
      tecplot_file << "\n";
    }
    
    /*--- Write the element connectivty in Tecplot format. ---*/
    for(iElem = 0; iElem < nElem[iZone]; iElem++) {
      if (elem[iZone][iElem]->GetVTK_Type() == TRIANGLE) {
        tecplot_file <<
        elem[iZone][iElem]->GetNode(0)+1 <<" "<< elem[iZone][iElem]->GetNode(1)+1 <<" "<<
        elem[iZone][iElem]->GetNode(2)+1 <<" "<< elem[iZone][iElem]->GetNode(2)+1 << endl;
      }
      if (elem[iZone][iElem]->GetVTK_Type() == RECTANGLE) {
        tecplot_file <<
        elem[iZone][iElem]->GetNode(0)+1 <<" "<< elem[iZone][iElem]->GetNode(1)+1 <<" "<<
        elem[iZone][iElem]->GetNode(2)+1 <<" "<< elem[iZone][iElem]->GetNode(3)+1 << endl;
      }
      if (elem[iZone][iElem]->GetVTK_Type() == TETRAHEDRON) {
        tecplot_file <<
        elem[iZone][iElem]->GetNode(0)+1 <<" "<< elem[iZone][iElem]->GetNode(1)+1 <<" "<<
        elem[iZone][iElem]->GetNode(2)+1 <<" "<< elem[iZone][iElem]->GetNode(2)+1 <<" "<<
        elem[iZone][iElem]->GetNode(3)+1 <<" "<< elem[iZone][iElem]->GetNode(3)+1 <<" "<<
        elem[iZone][iElem]->GetNode(3)+1 <<" "<< elem[iZone][iElem]->GetNode(3)+1 << endl;
      }
      if (elem[iZone][iElem]->GetVTK_Type() == HEXAHEDRON) {
        tecplot_file <<
        elem[iZone][iElem]->GetNode(0)+1 <<" "<< elem[iZone][iElem]->GetNode(1)+1 <<" "<<
        elem[iZone][iElem]->GetNode(2)+1 <<" "<< elem[iZone][iElem]->GetNode(3)+1 <<" "<<
        elem[iZone][iElem]->GetNode(4)+1 <<" "<< elem[iZone][iElem]->GetNode(5)+1 <<" "<<
        elem[iZone][iElem]->GetNode(6)+1 <<" "<< elem[iZone][iElem]->GetNode(7)+1 << endl;
      }
      if (elem[iZone][iElem]->GetVTK_Type() == PYRAMID) {
        tecplot_file <<
        elem[iZone][iElem]->GetNode(0)+1 <<" "<< elem[iZone][iElem]->GetNode(1)+1 <<" "<<
        elem[iZone][iElem]->GetNode(2)+1 <<" "<< elem[iZone][iElem]->GetNode(3)+1 <<" "<<
        elem[iZone][iElem]->GetNode(4)+1 <<" "<< elem[iZone][iElem]->GetNode(4)+1 <<" "<<
        elem[iZone][iElem]->GetNode(4)+1 <<" "<< elem[iZone][iElem]->GetNode(4)+1 << endl;
      }
      if (elem[iZone][iElem]->GetVTK_Type() == WEDGE) {
        tecplot_file <<
        elem[iZone][iElem]->GetNode(0)+1 <<" "<< elem[iZone][iElem]->GetNode(1)+1 <<" "<<
        elem[iZone][iElem]->GetNode(1)+1 <<" "<< elem[iZone][iElem]->GetNode(2)+1 <<" "<<
        elem[iZone][iElem]->GetNode(3)+1 <<" "<< elem[iZone][iElem]->GetNode(4)+1 <<" "<<
        elem[iZone][iElem]->GetNode(4)+1 <<" "<< elem[iZone][iElem]->GetNode(5)+1 << endl;
      }
    }
    
  }
  
  /*--- Close the Tecplot file after writing the final zone. ---*/
  tecplot_file.close();
  
}


void SetMeshFile(CPoint ***node, unsigned long *nPoint, CPrimalGrid ***elem,
                 unsigned long *nElem, CPrimalGrid**** bound, unsigned long **nElem_Bound,
                 CPrimalGrid**** newBound, unsigned long **nNewElem_Bound,
                 CConfig **config_container, unsigned short nZone, unsigned short nDim) {
  
  cout << "writing mesh files" << endl;
  
  ///////
  ////// Write the mesh files...
  //////
  
  unsigned long nPointNew, iElem_Bound, iElem, iPoint;
  unsigned short iNode, iZone, iDim, nMarker, iMarker, jMarker;
  ofstream su2_file;
  string Grid_Marker;
  
  char buffer[100];
  string filename;
  filename.assign(config_container[ZONE_0]->GetMesh_Out_FileName());
  strcpy (buffer, filename.c_str());
  
  
  /*--- Open .su2 grid file ---*/
  su2_file.precision(15);
  su2_file.open(buffer, ios::out);
  
  /*---  Write the total number of zones as the first line of the file. ---*/
  su2_file << "NZONE= " << nZone << endl;
  
  /*--- Loop over all zones. ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    
    /*--- Write the zone, dimension, and number of elements ---*/
    su2_file << "IZONE= " << iZone+1 << endl;
    su2_file << "NDIME= " << nDim << endl;
    su2_file << "NELEM= " << nElem[iZone] << endl;
    
    /*--- Write the connectivity for each element ---*/
    for (iElem = 0; iElem < nElem[iZone]; iElem++) {
      su2_file << elem[iZone][iElem]->GetVTK_Type();
      for (iNode = 0; iNode < elem[iZone][iElem]->GetnNodes(); iNode++)
        su2_file << "\t" << elem[iZone][iElem]->GetNode(iNode);
      su2_file << "\t"<<iElem<<endl;
    }
    
    /*--- Get the total number of markers in this zone and use the index
     to retrieve the number of new points that have been added to the
     zone (these are on the SEND boundary, which is indexed as nMarker+1). ---*/
    nMarker   = config_container[iZone]->GetnMarker_All();
    
    //NOTE that this is +1 bc we are using the old config, and the new recv boundary will be in that position
    nPointNew = nElem_Bound[iZone][nMarker+1];
    cout << nMarker <<" new points: "<<nPointNew<<endl;
    su2_file << "NPOIN= " << nPoint[iZone] << "\t" << nPoint[iZone] - nPointNew << endl;
    for (iPoint = 0; iPoint < nPoint[iZone]; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++)
        su2_file << scientific << "\t" << node[iZone][iPoint]->GetCoord(iDim) ;
      su2_file << "\t" << iPoint <<endl;
    }
    
    //////
    ////// WARNING: For now I am setting nMarker -= 2 so that the SEND/RECV markers are ignored by SU2_CFD.
    ////// Need to add interpolation...
    
    su2_file << "NMARK= " << nMarker+2 << endl;
    
    iMarker = 0;
    for (jMarker = 0; jMarker < config_container[iZone]->GetnMarker_All()+2; jMarker++) {
      
        if (bound[iZone][iMarker][0]->GetVTK_Type() != VERTEX) {
          
          
          Grid_Marker = config_container[iZone]->GetMarker_All_Tag(jMarker);
          su2_file << "MARKER_TAG= " << Grid_Marker <<endl;
          su2_file << "MARKER_ELEMS= " << nElem_Bound[iZone][iMarker] + nNewElem_Bound[iZone][jMarker] << endl;
          
          
          for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iZone][iMarker]; iElem_Bound++) {
            su2_file << bound[iZone][iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
            for (iNode = 0; iNode < bound[iZone][iMarker][iElem_Bound]->GetnNodes()-1; iNode++)
              su2_file << bound[iZone][iMarker][iElem_Bound]->GetNode(iNode) << "\t" ;
            iNode = bound[iZone][iMarker][iElem_Bound]->GetnNodes()-1;
            su2_file << bound[iZone][iMarker][iElem_Bound]->GetNode(iNode) << endl;
          }
          
          /*--- Write any new elements at the end of the list. ---*/
          if (nNewElem_Bound[iZone][jMarker] > 0) {
            for (iElem_Bound = 0; iElem_Bound < nNewElem_Bound[iZone][jMarker]; iElem_Bound++) {
              su2_file << newBound[iZone][jMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
              for (iNode = 0; iNode < newBound[iZone][jMarker][iElem_Bound]->GetnNodes()-1; iNode++)
                su2_file << newBound[iZone][jMarker][iElem_Bound]->GetNode(iNode) << "\t" ;
              iNode = newBound[iZone][jMarker][iElem_Bound]->GetnNodes()-1;
              su2_file << newBound[iZone][jMarker][iElem_Bound]->GetNode(iNode) << endl;
            }
          }
          
        }
        
        if (bound[iZone][iMarker][0]->GetVTK_Type() == VERTEX) {
          su2_file << "MARKER_TAG= SEND_RECEIVE" << endl;
          su2_file << "MARKER_ELEMS= " << nElem_Bound[iZone][iMarker]<< endl;
          if (config_container[iZone]->GetMarker_All_SendRecv(iMarker) > 0)
            su2_file << "SEND_TO= " << config_container[iZone]->GetMarker_All_SendRecv(iMarker) << endl;
          if (config_container[iZone]->GetMarker_All_SendRecv(iMarker) < 0)
            su2_file << "SEND_TO= " << config_container[iZone]->GetMarker_All_SendRecv(iMarker) << endl;
          
          //          su2_file << "DOMAIN= " << config->GetSlideDonor_Zone(config->GetMarker_All_Tag(iMarker)) << endl;
          //          su2_file << "SLIDING_MARKER= " << config->GetMarker_Sliding_Donor(config->GetMarker_All_Tag(iMarker)) << endl;
          
          for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iZone][iMarker]; iElem_Bound++) {
            su2_file << bound[iZone][iMarker][iElem_Bound]->GetVTK_Type() << "\t" <<
            bound[iZone][iMarker][iElem_Bound]->GetNode(0) << "\t" <<
            bound[iZone][iMarker][iElem_Bound]->GetRotation_Type() << "\t" <<
            bound[iZone][iMarker][iElem_Bound]->GetMatching_Zone() << endl;
          }
        }
        iMarker++;
      }
    
    su2_file << "NPERIODIC= " << 1 << endl;
    
    /*--- Periodic 0 correspond with no movement of the surface ---*/
    su2_file << "PERIODIC_INDEX= 0" << endl;
    su2_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << endl;
    su2_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << endl;
    su2_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << endl;
    su2_file << "\n" << endl;
    
  }
  
  su2_file.close();
  
}
