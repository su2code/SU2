/*!
 * \file geometry_structure.hpp
 * \brief Headers of the main subroutines for creating the geometrical structure.
 *        The subroutines and functions are in the <i>geometry_structure.cpp</i> file.
 * \author F. Palacios, T. Economon
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

#pragma once

#include "./mpi_structure.hpp"

#ifdef HAVE_METIS
  #include "metis.h"
#endif
#ifdef HAVE_PARMETIS
extern "C" {
#include "parmetis.h"
}
#endif
#ifdef HAVE_CGNS
  #include "fem_cgns_elements.hpp"
#endif
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <climits>

#include "primal_grid_structure.hpp"
#include "dual_grid_structure.hpp"
#include "config_structure.hpp"
#include "fem_standard_element.hpp"

using namespace std;

/*!
 * \class CUnsignedLong2T
 * \brief Help class used to store two unsigned longs as one entity.
 */
class CUnsignedLong2T {
public:
  unsigned long long0;  /*!< \brief First long to store in this class. */
  unsigned long long1;  /*!< \brief Second long to store in this class. */

  /* Constructors and destructors. */
  CUnsignedLong2T();
  ~CUnsignedLong2T();

  CUnsignedLong2T(const unsigned long a, const unsigned long b);

  CUnsignedLong2T(const CUnsignedLong2T &other);

  /* Operators. */
  CUnsignedLong2T& operator=(const CUnsignedLong2T &other);

  bool operator<(const CUnsignedLong2T &other) const;

  bool operator==(const CUnsignedLong2T &other) const;

private:
  /* Copy function. */
  void Copy(const CUnsignedLong2T &other);
};

/*!
 * \class CUnsignedShort2T
 * \brief Help class used to store two unsigned shorts as one entity.
 */
class CUnsignedShort2T {
public:
  unsigned short short0;  /*!< \brief First short to store in this class. */
  unsigned short short1;  /*!< \brief Second short to store in this class. */

  /* Constructors and destructors. */
  CUnsignedShort2T();
  ~CUnsignedShort2T();

  CUnsignedShort2T(const unsigned short a, const unsigned short b);

  CUnsignedShort2T(const CUnsignedShort2T &other);

  /* Operators. */
  CUnsignedShort2T& operator=(const CUnsignedShort2T &other);

  bool operator<(const CUnsignedShort2T &other) const;

  bool operator==(const CUnsignedShort2T &other) const;

private:
  /* Copy function. */
  void Copy(const CUnsignedShort2T &other);
};

/*!
 * \class CFaceOfElement
 * \brief Class used in the partitioning of the FEM grid as well as the building of
          the faces of DG. It stores a face of an element.
 */
class CFaceOfElement {
public:
  unsigned short nCornerPoints;          /*!< \brief Number of corner points of the face. */
  unsigned long  cornerPoints[4];        /*!< \brief Global ID's of ther corner points. */
  unsigned long  elemID0, elemID1;       /*!< \brief Element ID's to the left and right. */
  unsigned short nPolyGrid0, nPolyGrid1; /*!< \brief Polynomial degrees of the grid of the elements
                                                     to the left and right. */
  unsigned short nPolySol0,  nPolySol1;  /*!< \brief Polynomial degrees of the solution of the elements
                                                     to the left and right. */
  unsigned short nDOFsElem0, nDOFsElem1; /*!< \brief Number of DOFs of the elements to the left and right. */
  unsigned short elemType0,  elemType1;  /*!< \brief Type of the elements to the left and right. */
  unsigned short faceID0, faceID1;       /*!< \brief The local face ID in the corresponding elements
                                                     to the left and right of the face. */
  unsigned short periodicIndex;          /*!< \brief Periodic indicator of the face. A value of 0 means no
                                                     periodic face. A value larger than 0 gives the index of
                                                     the periodic boundary + 1. */
  unsigned short periodicIndexDonor;     /*!< \brief Periodic indicator of the donor face. A value of 0 means no
                                                     periodic donor face. A value larger than 0 gives the index of
                                                     the periodic donor boundary + 1. */
  short faceIndicator;                   /*!< \brief The corresponding boundary marker if this face is on a
                                                     boundary. For an internal face the value is -1,
                                                     while an invalidated face has the value -2. */
  bool JacFaceIsConsideredConstant;      /*!< \brief Whether or not the Jacobian of the transformation
                                                     to the standard element is considered constant. */
  bool elem0IsOwner;                     /*!< \brief Whether or not the neighboring element 0 is the owner
                                                     of the face. If false, element 1 is the owner. */

  /* Standard constructor and destructor. */
  CFaceOfElement();
  ~CFaceOfElement(){}

  /* Alternative constructor to set the corner points. */
  CFaceOfElement(const unsigned short VTK_Type,
                 const unsigned short nPoly,
                 const unsigned long  *Nodes);

  /* Copy constructor and assignment operator. */
  CFaceOfElement(const CFaceOfElement &other);

  CFaceOfElement& operator=(const CFaceOfElement &other);

  /* Less than operator. Needed for the sorting and searching. */
  bool operator<(const CFaceOfElement &other) const;

  /* Equal operator. Needed for removing double entities. */
  bool operator ==(const CFaceOfElement &other) const;

  /*--- Member function, which creates a unique numbering for the corner points.
        A sort in increasing order is OK for this purpose.                       ---*/
  void CreateUniqueNumbering(void);

  /*--- Member function, which creates a unique numbering for the corner points
        while the orientation is taken into account. ---*/
  void CreateUniqueNumberingWithOrientation(void);

private:
  /*--- Copy function, which copies the data of the given object into the current object. ---*/
  void Copy(const CFaceOfElement &other);
};

/*!
 * \class CBoundaryFace
 * \brief Help class used in the partitioning of the FEM grid.
          It stores a boundary element.
 */
class CBoundaryFace {
 public:
  unsigned short VTK_Type, nPolyGrid, nDOFsGrid;
  unsigned long  globalBoundElemID, domainElementID;
  vector<unsigned long>  Nodes;

  /* Standard constructor and destructor. Nothing to be done. */
  CBoundaryFace(){}
  ~CBoundaryFace(){}

  /* Copy constructor and assignment operator. */
  CBoundaryFace(const CBoundaryFace &other);

  CBoundaryFace& operator=(const CBoundaryFace &other);

  /* Less than operator. Needed for the sorting. */
  bool operator<(const CBoundaryFace &other) const;

private:
  /*--- Copy function, which copies the data of the given object into the current object. ---*/
  void Copy(const CBoundaryFace &other);
};

/*!
 * \class CMatchingFace
 * \brief Help class used to determine whether or not (periodic) faces match.
 */
class CMatchingFace {
public:
  unsigned short nCornerPoints;          /*!< \brief Number of corner points of the face. */
  unsigned short nDim;                   /*!< \brief Number of spatial dimensions. */
  unsigned short nPoly;                  /*!< \brief Polynomial degree of the face. */
  unsigned short nDOFsElem;              /*!< \brief Number of DOFs of the relevant adjacent element. */
  unsigned short elemType;               /*!< \brief Type of the adjacent element. */
  unsigned long  elemID;                 /*!< \brief The relevant adjacent element ID. */
  su2double cornerCoor[4][3];            /*!< \brief Coordinates of the corner points of the face. */
  su2double tolForMatching;              /*!< \brief Tolerance for this face for matching points. */

  /* Standard constructor. */
  CMatchingFace();

  /* Destructor, nothing to be done. */
  ~CMatchingFace(){}

  /* Copy constructor and assignment operator. */
  CMatchingFace(const CMatchingFace &other);

  CMatchingFace& operator=(const CMatchingFace &other);

  /* Less than operator. Needed for the sorting and searching. */
  bool operator<(const CMatchingFace &other) const;

  /*--- Member function, which sorts the coordinates of the face. ---*/
  void SortFaceCoordinates(void);

private:
  /*--- Copy function, which copies the data of the given object into the current object. ---*/
  void Copy(const CMatchingFace &other);
};

/*!
 * \class CGeometry
 * \brief Parent class for defining the geometry of the problem (complete geometry, 
 *        multigrid agglomerated geometry, only boundary geometry, etc..)
 * \author F. Palacios
 */
class CGeometry {
protected:
  int rank, 	/*!< \brief MPI Rank. */
  size;       	/*!< \brief MPI Size. */
	unsigned long nPoint,	/*!< \brief Number of points of the mesh. */
	nPointDomain,						/*!< \brief Number of real points of the mesh. */
	nPointGhost,					/*!< \brief Number of ghost points of the mesh. */
	nPointNode,					/*!< \brief Size of the node array allocated to hold CPoint objects. */
  Global_nPoint,	/*!< \brief Total number of nodes in a simulation across all processors (including halos). */
	Global_nPointDomain,	/*!< \brief Total number of nodes in a simulation across all processors (excluding halos). */
	nElem,					/*!< \brief Number of elements of the mesh. */
  Global_nElem,	/*!< \brief Total number of elements in a simulation across all processors (all types). */
  Global_nElemDomain,  /*!< \brief Total number of elements in a simulation across all processors (excluding halos). */
	nEdge,					/*!< \brief Number of edges of the mesh. */
	nFace,					/*!< \brief Number of faces of the mesh. */
  nelem_edge,             /*!< \brief Number of edges in the mesh. */
  Global_nelem_edge,      /*!< \brief Total number of edges in the mesh across all processors. */
  nelem_triangle,       /*!< \brief Number of triangles in the mesh. */
  Global_nelem_triangle,       /*!< \brief Total number of triangles in the mesh across all processors. */
  nelem_quad,           /*!< \brief Number of quadrangles in the mesh. */
  Global_nelem_quad,           /*!< \brief Total number of quadrangles in the mesh across all processors. */
  nelem_tetra,          /*!< \brief Number of tetrahedra in the mesh. */
  Global_nelem_tetra,          /*!< \brief Total number of tetrahedra in the mesh across all processors. */
  nelem_hexa,           /*!< \brief Number of hexahedra in the mesh. */
  Global_nelem_hexa,           /*!< \brief Total number of hexahedra in the mesh across all processors. */
  nelem_prism,          /*!< \brief Number of prisms in the mesh. */
  Global_nelem_prism,          /*!< \brief Total number of prisms in the mesh across all processors. */
  nelem_pyramid,        /*!< \brief Number of pyramids in the mesh. */
  Global_nelem_pyramid,        /*!< \brief Total number of pyramids in the mesh across all processors. */
  nelem_edge_bound,           /*!< \brief Number of edges on the mesh boundaries. */
  Global_nelem_edge_bound,           /*!< \brief Total number of edges on the mesh boundaries across all processors. */
  nelem_triangle_bound,          /*!< \brief Number of triangles on the mesh boundaries. */
  Global_nelem_triangle_bound,          /*!< \brief Total number of triangles on the mesh boundaries across all processors. */
  nelem_quad_bound,        /*!< \brief Number of quads on the mesh boundaries. */
  Global_nelem_quad_bound;        /*!< \brief Total number of quads on the mesh boundaries across all processors. */
	unsigned short nDim,	/*!< \brief Number of dimension of the problem. */
	nZone,								/*!< \brief Number of zones in the problem. */
	nMarker;				/*!< \brief Number of different markers of the mesh. */
  unsigned short MGLevel;         /*!< \brief The mesh level index for the current geometry container. */
  unsigned long Max_GlobalPoint;  /*!< \brief Greater global point in the domain local structure. */

  /* --- Custom boundary variables --- */
  su2double **CustomBoundaryTemperature;
  su2double **CustomBoundaryHeatFlux;

public:
	unsigned long *nElem_Bound;			/*!< \brief Number of elements of the boundary. */
	string *Tag_to_Marker;	/*!< \brief If you know the index of the boundary (depend of the 
							 grid definition), it gives you the maker (where the boundary 
							 is stored from 0 to boundaries). */	
	CPrimalGrid** elem;	/*!< \brief Element vector (primal grid information). */
	CPrimalGrid** face;			/*!< \brief Face vector (primal grid information). */
	CPrimalGrid*** bound;	/*!< \brief Boundary vector (primal grid information). */
	CPoint** node;			/*!< \brief Node vector (dual grid information). */
	CEdge** edge;			/*!< \brief Edge vector (dual grid information). */
	CVertex*** vertex;		/*!< \brief Boundary Vertex vector (dual grid information). */
  CTurboVertex**** turbovertex; /*!< \brief Boundary Vertex vector ordered for turbomachinery calculation(dual grid information). */
  unsigned long *nVertex;	/*!< \brief Number of vertex for each marker. */
  unsigned short *nSpanWiseSections; /*!< \brief Number of Span wise section for each turbo marker, indexed by inflow/outflow */
  unsigned short *nSpanSectionsByMarker; /*! <\brief Number of Span wise section for each turbo marker, indexed by marker.  Needed for deallocation.*/
  unsigned short nTurboPerf; /*!< \brief Number of Span wise section for each turbo marker. */
  su2double **SpanWiseValue; /*!< \brief Span wise values for each turbo marker. */
  long **nVertexSpan; /*! <\brief number of vertexes for span wise section for each marker.  */
  unsigned long **nTotVertexSpan; /*! <\brief number of vertexes at each span wise section for each marker.  */
  unsigned long nVertexSpanMax[3]; /*! <\brief max number of vertexes for each span section for each marker flag.  */
  su2double ***AverageTurboNormal; /*! <\brief Average boundary normal at each span wise section for each marker in the turbomachinery frame of reference.*/
  su2double ***AverageNormal; /*! <\brief Average boundary normal at each span wise section for each marker.*/
  su2double ***AverageGridVel; /*! <\brief Average boundary grid velocity at each span wise section for each marker.*/
  su2double **AverageTangGridVel; /*! <\brief Average tangential rotational speed at each span wise section for each marker.*/
  su2double **SpanArea; /*! <\brief Area at each span wise section for each marker.*/
  su2double **MaxAngularCoord; /*! <\brief Max angular pitch at each span wise section for each marker.*/
  su2double **MinAngularCoord; /*! <\brief Max angular pitch at each span wise section for each marker.*/
  su2double **MinRelAngularCoord; /*! <\brief Min relative angular coord at each span wise section for each marker.*/
  su2double **TurboRadius; /*! <\brief Radius at each span wise section for each marker.*/
  su2double **TangGridVelIn, **TangGridVelOut; /*! <\brief Average tangential rotational speed at each span wise section for each turbomachinery marker.*/
  su2double **SpanAreaIn, **SpanAreaOut; /*! <\brief Area at each span wise section for each turbomachinery marker.*/
  su2double **TurboRadiusIn, **TurboRadiusOut; /*! <\brief Radius at each span wise section for each turbomachinery marker*/

  unsigned short nCommLevel;		/*!< \brief Number of non-blocking communication levels. */
	vector<unsigned long> PeriodicPoint[MAX_NUMBER_PERIODIC][2];			/*!< \brief PeriodicPoint[Periodic bc] and return the point that
																			 must be sent [0], and the image point in the periodic bc[1]. */
	vector<unsigned long> PeriodicElem[MAX_NUMBER_PERIODIC];				/*!< \brief PeriodicElem[Periodic bc] and return the elements that 
																			 must be sent. */
  
  short *Marker_All_SendRecv;
  
	/*--- Create vectors and distribute the values among the different planes queues ---*/
	vector<vector<su2double> > Xcoord_plane; /*!< \brief Vector containing x coordinates of new points appearing on a single plane */
	vector<vector<su2double> > Ycoord_plane; /*!< \brief Vector containing y coordinates of  new points appearing on a single plane */
	vector<vector<su2double> > Zcoord_plane; 	/*!< \brief Vector containing z coordinates of  new points appearing on a single plane */
	vector<vector<su2double> > FaceArea_plane; /*!< \brief Vector containing area/volume associated with  new points appearing on a single plane */
	vector<vector<unsigned long> > Plane_points; /*!< \brief Vector containing points appearing on a single plane */

	vector<su2double> XCoordList;	/*!< \brief Vector containing points appearing on a single plane */
	CPrimalGrid*** newBound;            /*!< \brief Boundary vector for new periodic elements (primal grid information). */
	unsigned long *nNewElem_Bound;			/*!< \brief Number of new periodic elements of the boundary. */
  
  /*--- Partitioning-specific variables ---*/
  map<unsigned long,unsigned long> Global_to_Local_Elem;
  unsigned long xadj_size;
  unsigned long adjacency_size;
  unsigned long *starting_node;
  unsigned long *ending_node;
  unsigned long *npoint_procs;
  unsigned long *nPoint_Linear;
#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
  idx_t * adjacency;
  idx_t * xadj;
#endif
#endif
  
  /*--- Data structures for point-to-point MPI communications. ---*/

  int countPerPoint;                  /*!< \brief Maximum number of pieces of data sent per vertex in point-to-point comms. */
  int nP2PSend;                       /*!< \brief Number of sends during point-to-point comms. */
  int nP2PRecv;                       /*!< \brief Number of receives during point-to-point comms. */
  int *nPoint_P2PSend;                /*!< \brief Data structure holding number of vertices for each send in point-to-point comms. */
  int *nPoint_P2PRecv;                /*!< \brief Data structure holding number of vertices for each recv in point-to-point comms. */
  int *Neighbors_P2PSend;             /*!< \brief Data structure holding the ranks of the neighbors for point-to-point send comms. */
  int *Neighbors_P2PRecv;             /*!< \brief Data structure holding the ranks of the neighbors for point-to-point recv comms. */
  map<int, int> P2PSend2Neighbor;     /*!< \brief Data structure holding the reverse mapping of the ranks of the neighbors for point-to-point send comms. */
  map<int, int> P2PRecv2Neighbor;     /*!< \brief Data structure holding the reverse mapping of the ranks of the neighbors for point-to-point recv comms. */
  unsigned long *Local_Point_P2PSend; /*!< \brief Data structure holding the local index of all vertices to be sent in point-to-point comms. */
  unsigned long *Local_Point_P2PRecv; /*!< \brief Data structure holding the local index of all vertices to be received in point-to-point comms. */
  su2double *bufD_P2PRecv;            /*!< \brief Data structure for su2double point-to-point receive. */
  su2double *bufD_P2PSend;            /*!< \brief Data structure for su2double point-to-point send. */
  unsigned short *bufS_P2PRecv;       /*!< \brief Data structure for unsigned long point-to-point receive. */
  unsigned short *bufS_P2PSend;       /*!< \brief Data structure for unsigned long point-to-point send. */
  SU2_MPI::Request *req_P2PSend;      /*!< \brief Data structure for point-to-point send requests. */
  SU2_MPI::Request *req_P2PRecv;      /*!< \brief Data structure for point-to-point recv requests. */

  /*--- Data structures for periodic communications. ---*/

  int countPerPeriodicPoint;                /*!< \brief Maximum number of pieces of data sent per vertex in periodic comms. */
  int nPeriodicSend;                        /*!< \brief Number of sends during periodic comms. */
  int nPeriodicRecv;                        /*!< \brief Number of receives during periodic comms. */
  int *nPoint_PeriodicSend;                 /*!< \brief Data structure holding number of vertices for each send in periodic comms. */
  int *nPoint_PeriodicRecv;                 /*!< \brief Data structure holding number of vertices for each recv in periodic comms. */
  int *Neighbors_PeriodicSend;              /*!< \brief Data structure holding the ranks of the neighbors for periodic send comms. */
  int *Neighbors_PeriodicRecv;              /*!< \brief Data structure holding the ranks of the neighbors for periodic recv comms. */
  map<int, int> PeriodicSend2Neighbor;      /*!< \brief Data structure holding the reverse mapping of the ranks of the neighbors for periodic send comms. */
  map<int, int> PeriodicRecv2Neighbor;      /*!< \brief Data structure holding the reverse mapping of the ranks of the neighbors for periodic recv comms. */
  unsigned long *Local_Point_PeriodicSend;  /*!< \brief Data structure holding the local index of all vertices to be sent in periodic comms. */
  unsigned long *Local_Point_PeriodicRecv;  /*!< \brief Data structure holding the local index of all vertices to be received in periodic comms. */
  unsigned long *Local_Marker_PeriodicSend; /*!< \brief Data structure holding the local index of the periodic marker for a particular vertex to be sent in periodic comms. */
  unsigned long *Local_Marker_PeriodicRecv; /*!< \brief Data structure holding the local index of the periodic marker for a particular vertex to be received in periodic comms. */
  su2double *bufD_PeriodicRecv;             /*!< \brief Data structure for su2double periodic receive. */
  su2double *bufD_PeriodicSend;             /*!< \brief Data structure for su2double periodic send. */
  unsigned short *bufS_PeriodicRecv;        /*!< \brief Data structure for unsigned long periodic receive. */
  unsigned short *bufS_PeriodicSend;        /*!< \brief Data structure for unsigned long periodic send. */
  SU2_MPI::Request *req_PeriodicSend;       /*!< \brief Data structure for periodic send requests. */
  SU2_MPI::Request *req_PeriodicRecv;       /*!< \brief Data structure for periodic recv requests. */
  
	/*!
	 * \brief Constructor of the class.
	 */
	CGeometry(void);

	/*! 
	 * \brief Destructor of the class.
	 */
	virtual ~CGeometry(void);
  
  /*!
   * \brief Routine to launch non-blocking recvs only for all periodic communications. Note that this routine is called by any class that has loaded data into the generic communication buffers.
   * \brief Routine to set up persistent data structures for point-to-point MPI communications.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void PreprocessP2PComms(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Routine to allocate buffers for point-to-point MPI communications. Also called to dynamically reallocate if not enough memory is found for comms during runtime.
   * \param[in] val_countPerPoint - Maximum count of the data type per vertex in point-to-point comms, e.g., nPrimvarGrad*nDim.
   */
  void AllocateP2PComms(unsigned short val_countPerPoint);
  
  /*!
   * \brief Routine to launch non-blocking recvs only for all point-to-point communication with neighboring partitions. Note that this routine is called by any class that has loaded data into the generic communication buffers.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be communicated.
   * \param[in] val_reverse  - Boolean controlling forward or reverse communication between neighbors.
   */
  void PostP2PRecvs(CGeometry *geometry, CConfig *config, unsigned short commType, bool val_reverse);
  
  /*!
   * \brief Routine to launch a single non-blocking send once the buffer is loaded for a point-to-point commucation. Note that this routine is called by any class that has loaded data into the generic communication buffers.
   * \param[in] geometry     - Geometrical definition of the problem.
   * \param[in] config       - Definition of the particular problem.
   * \param[in] commType     - Enumerated type for the quantity to be communicated.
   * \param[in] val_iMessage - Index of the message in the order they are stored.
   * \param[in] val_reverse  - Boolean controlling forward or reverse communication between neighbors.
   */
  void PostP2PSends(CGeometry *geometry, CConfig *config, unsigned short commType, int val_iMessage, bool val_reverse);
  
  /*!
   * \brief Routine to set up persistent data structures for periodic communications.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void PreprocessPeriodicComms(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Routine to allocate buffers for periodic communications. Also called to dynamically reallocate if not enough memory is found for comms during runtime.
   * \param[in] val_countPerPeriodicPoint - Maximum count of the data type per vertex in periodic comms, e.g., nPrimvarGrad*nDim.
   */
  void AllocatePeriodicComms(unsigned short val_countPerPeriodicPoint);

  /*!
   * \brief Routine to launch non-blocking recvs only for all periodic communication with neighboring partitions. Note that this routine is called by any class that has loaded data into the generic communication buffers.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be communicated.
   */
  void PostPeriodicRecvs(CGeometry *geometry, CConfig *config, unsigned short commType);
  
  /*!
   * \brief Routine to launch a single non-blocking send once the buffer is loaded for a periodic commucation. Note that this routine is called by any class that has loaded data into the generic communication buffers.
   * \param[in] geometry     - Geometrical definition of the problem.
   * \param[in] config       - Definition of the particular problem.
   * \param[in] commType     - Enumerated type for the quantity to be communicated.
   * \param[in] val_iMessage - Index of the message in the order they are stored.
   */
  void PostPeriodicSends(CGeometry *geometry, CConfig *config, unsigned short commType, int val_iMessage);

  /*!
   * \brief Routine to load a geometric quantity into the data structures for MPI point-to-point communication and to launch non-blocking sends and recvs for all point-to-point communication with neighboring partitions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be communicated.
   */
  void InitiateComms(CGeometry *geometry, CConfig *config, unsigned short commType);
  
  /*!
   * \brief Routine to complete the set of non-blocking communications launched by InitiateComms() and unpacking of the data into the geometry class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be unpacked.
   */
  void CompleteComms(CGeometry *geometry, CConfig *config, unsigned short commType);
  
	/*! 
	 * \brief Get number of coordinates.
	 * \return Number of coordinates.
	 */
	unsigned short GetnDim(void);

	/*! 
	 * \brief Get number of zones.
	 * \return Number of zones.
	 */
	unsigned short GetnZone(void);

	/*! 
	 * \brief Get number of points.
	 * \return Number of points.
	 */
	unsigned long GetnPoint(void);

	/*! 
	 * \brief Get number of real points (that belong to the domain).
	 * \return Number of real points.
	 */
	unsigned long GetnPointDomain(void);

  /*!
	 * \brief Get number of elements.
	 * \return Number of elements.
	 */
	unsigned long GetnLine(void);
  
	/*! 
	 * \brief Get number of elements.
	 * \return Number of elements.
	 */
	unsigned long GetnElem(void);
  
	/*! 
	 * \brief Get number of edges.
	 * \return Number of edges.
	 */
	unsigned long GetnEdge(void);

	/*! 
	 * \brief Get number of markers.
	 * \return Number of markers.
	 */
	unsigned short GetnMarker(void);

	/*!
	 * \brief Get number of vertices.
	 * \param[in] val_marker - Marker of the boundary.
	 * \return Number of vertices.
	 */
	su2double* GetSpanWiseValue(unsigned short val_marker);

	/*! 
	 * \brief Get number of vertices.
	 * \param[in] val_marker - Marker of the boundary.
	 * \return Number of vertices.
	 */
	unsigned long GetnVertex(unsigned short val_marker);

	/*!
	 * \brief Get number of span wise section.
	 * \param[in] marker_flag - flag of the turbomachinery boundary.
	 * \return Number of span wise section.
	 */
	unsigned short GetnSpanWiseSections(unsigned short marker_flag);

	/*!
	 * \brief Get number of vertices.
	 * \param[in] val_marker - Marker of the boundary.
	 * \return Number of vertices.
	 */
	unsigned long GetnVertexSpan(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Get number of frequencies per span for NRBC.
	 * \param[in] val_marker - Marker of the boundary.
	 * \return Number of frequencies for NRBC.
	 */
	unsigned long GetnFreqSpan(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Get number of vertices.
	 * \param[in] val_marker - Marker of the boundary.
	 * \return Number of vertices.
	 */
	unsigned long GetnVertexSpanMax(unsigned short marker_flag);

	/*!
	 * \brief Get number of max frequencies for initializing the Fourier Coefficient for NR BC.
	 * \param[in] marker_flag - Marker of the boundary.
	 * \return Number of frequencies.
	 */
	unsigned long GetnFreqSpanMax(unsigned short marker_flag);

	/*!
	 * \brief Get number of vertices.
	 * \param[in] val_marker - Marker of the boundary.
	 * \return Number of vertices.
	 */
	void SetnVertexSpanMax(unsigned short marker_flag, unsigned long nVertMax);

	/*! 
	 * \brief Get the edge index from using the nodes of the edge.
	 * \param[in] first_point - First point of the edge.
	 * \param[in] second_point - Second point of the edge.
	 * \return Index of the edge.
	 */		
	long FindEdge(unsigned long first_point, unsigned long second_point);

    /*!
	 * \brief Get the edge index from using the nodes of the edge.
	 * \param[in] first_point - First point of the edge.
	 * \param[in] second_point - Second point of the edge.
	 * \return Index of the edge.
	 */
	bool CheckEdge(unsigned long first_point, unsigned long second_point);
    
	/*! 
	 * \brief Get the distance between a plane (defined by three point) and a point.
	 * \param[in] Coord - Coordinates of the point.
	 * \param[in] iCoord - Coordinates of the first point that defines the plane.
	 * \param[in] jCoord - Coordinates of the second point that defines the plane.
	 * \param[in] kCoord - Coordinates of the third point that defines the plane.
	 * \return Signed distance.
	 */		
	su2double Point2Plane_Distance(su2double *Coord, su2double *iCoord, su2double *jCoord, su2double *kCoord);

	/*! 
	 * \brief Create a file for testing the geometry.
	 */		
	void TestGeometry(void);

	/*! 
	 * \brief A virtual member.
	 * \param[in] val_nmarker - Number of markers.
	 */
	void SetnMarker(unsigned short val_nmarker);

	/*! 
	 * \brief Set the number of dimensions of the problem.
	 * \param[in] val_nDim - Number of dimensions.
	 */
	void SetnDim(unsigned short val_nDim);

	/*! 
	 * \brief Get the index of a marker.
	 * \param[in] val_marker - Marker of the boundary.
	 * \return Index of the marker in the grid defintion.
	 */	
	string GetMarker_Tag(unsigned short val_marker);

	/*! 
	 * \brief Set index of a marker.
	 * \param[in] val_marker - Marker of the boundary.
	 * \param[in] val_index - Index of the marker.
	 */		
	void SetMarker_Tag(unsigned short val_marker, string val_index);

	/*! 
	 * \brief Set the number of boundary elements.
	 * \param[in] val_marker - Marker of the boundary.
	 * \param[in] val_nelem_bound - Number of boundary elements.
	 */	
	void SetnElem_Bound(unsigned short val_marker, unsigned long val_nelem_bound);

	/*! 
	 * \brief Set the number of grid points.
	 * \param[in] val_npoint - Number of grid points.
	 */	
	void SetnPoint(unsigned long val_npoint);

	/*! 
	 * \brief Set the number of grid points in the domain.
	 * \param[in] val_npoint - Number of grid points in the domain.
	 */	
	void SetnPointDomain(unsigned long val_npoint);

	/*! 
	 * \brief Set the number of grid elements.
	 * \param[in] val_nelem - Number of grid elements.
	 */
	void SetnElem(unsigned long val_nelem);

	/*! 
	 * \brief Get the number of boundary elements.
	 * \param[in] val_marker - Marker of the boundary.
	 */
	unsigned long GetnElem_Bound(unsigned short val_marker);

  /*!
	 * \brief Get the number of elements in vtk fortmat.
	 */
	unsigned long GetMax_GlobalPoint(void);

	/*! 
	 * \brief A virtual function.
	 * \param[in] first_elem - Identification of the first element.
	 * \param[in] second_elem - Identification of the second element.
	 * \param[in] face_first_elem - Index of the common face for the first element.
	 * \param[in] face_second_elem - Index of the common face for the second element.
	 */
	virtual bool FindFace(unsigned long first_elem, unsigned long second_elem, unsigned short &face_first_elem, 
			unsigned short &face_second_elem);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.		 
	 */
	virtual void ComputeWall_Distance(CConfig *config);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.		 
	 */
	virtual void SetPositive_ZArea(CConfig *config);
  
	/*! 
	 * \brief A virtual member.
	 */	
	virtual void SetPoint_Connectivity(void);
  
  /*!
	 * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetRCM_Ordering(CConfig *config);
  
	/*!
	 * \brief A virtual member.
	 */		
	virtual void SetElement_Connectivity(void);

	/*! 
	 * \brief A virtual member.
	 */
	void SetEdges(void);

	/*! 
	 * \brief A virtual member.
	 */
	virtual void SetFaces(void);

	/*!
	 * \brief A virtual member.
	 */
	virtual void SetBoundVolume(void);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetVertex(CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeNSpan(CConfig *config, unsigned short val_iZone, unsigned short marker_flag, bool allocate);

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetTurboVertex(CConfig *config, unsigned short val_iZone, unsigned short marker_flag, bool allocate);

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void UpdateTurboVertex(CConfig *config, unsigned short val_iZone, unsigned short marker_flag);

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetAvgTurboValue(CConfig *config, unsigned short val_iZone, unsigned short marker_flag, bool allocate);

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void GatherInOutAverageValues(CConfig *config, bool allocate);

  /*!
   * \brief A virtual member.
   */
  virtual void SetVertex(void);

  /*!
   * \brief A virtual member.
   */
  virtual void SetCoord_CG(void);

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetMaxLength(CConfig* config);

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   * \param[in] action - Allocate or not the new elements.
   */
  virtual void SetControlVolume(CConfig *config, unsigned short action);

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   * \param[in] action - Allocate or not the new elements.
   */
  virtual void VisualizeControlVolume(CConfig *config, unsigned short action);
  
	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void MatchNearField(CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void MatchActuator_Disk(CConfig *config);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void MatchInterface(CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void MatchPeriodic(CConfig *config, unsigned short val_periodic);
  
	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] action - Allocate or not the new elements.		 
	 */
	virtual void SetBoundControlVolume(CConfig *config, unsigned short action);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config_filename - Name of the file where the tecplot information is going to be stored.
	 */
	virtual void SetTecPlot(char config_filename[MAX_STRING_SIZE], bool new_file);

	/*! 
	 * \brief A virtual member.
   * \param[in] mesh_filename - Name of the file where the tecplot information is going to be stored.
   * \param[in] new_file - Boolean to decide if aopen a new file or add to a old one
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetBoundTecPlot(char mesh_filename[MAX_STRING_SIZE], bool new_file, CConfig *config);
  
	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.		 
	 */
	virtual void Check_IntElem_Orientation(CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Check_BoundElem_Orientation(CConfig *config);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.		 
	 */
	virtual void SetColorGrid(CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetColorGrid_Parallel(CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetColorFEMGrid_Parallel(CConfig *config);

  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
  virtual void DivideConnectivity(CConfig *config, unsigned short Elem_Type);

	/*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_domain - Number of domains for parallelization purposes.		 
	 */
	virtual void SetSendReceive(CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_domain - Number of domains for parallelization purposes.
	 */
	virtual void SetBoundaries(CConfig *config);
  
	/*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */	
	virtual void SetCoord(CGeometry *geometry);

        /*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
         * \param[in] val_marker - Index of the boundary marker.
	 */
        virtual void SetMultiGridWallHeatFlux(CGeometry *geometry, unsigned short val_marker);

        /*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
         * \param[in] val_marker - Index of the boundary marker.
	 */
        virtual void SetMultiGridWallTemperature(CGeometry *geometry, unsigned short val_marker);

	/*! 
	 * \brief A virtual member.
	 * \param[in] val_nSmooth - Number of smoothing iterations.
	 * \param[in] val_smooth_coeff - Relaxation factor.
	 * \param[in] config - Definition of the particular problem.
	 */	
	virtual void SetCoord_Smoothing(unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config);

	/*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */	
	virtual void SetPoint_Connectivity(CGeometry *geometry);

	/*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */	
	virtual void SetVertex(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] action - Allocate or not the new elements.		 
	 */	
	virtual void SetControlVolume(CConfig *config, CGeometry *geometry, unsigned short action);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] action - Allocate or not the new elements.		 
	 */	
	virtual void SetBoundControlVolume(CConfig *config, CGeometry *geometry, unsigned short action);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_mesh_out_filename - Name of the output file.
	 */	
	virtual void SetMeshFile(CConfig *config, string val_mesh_out_filename);
  
    /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_mesh_out_filename - Name of the output file.
	 */
	virtual void SetMeshFile(CGeometry *geometry, CConfig *config, string val_mesh_out_filename);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetBoundSensitivity(CConfig *config);

  /*!
   * \brief Set the data containers for customized boundary conditions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetCustomBoundary(CConfig *config);


  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iZone - Index of the current zone.
   */
  virtual void SetRotationalVelocity(CConfig *config, unsigned short val_iZone, bool print);

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetShroudVelocity(CConfig *config);

   /*!
    * \brief A virtual member.
    * \param[in] config - Definition of the particular problem.
    */
   virtual void SetTranslationalVelocity(CConfig *config, unsigned short val_iZone, bool print);

   /*!
    * \brief A virtual member.
    * \param[in] config - Definition of the particular problem.
    * \param[in] iter - Current physical time step.
    */
   virtual void SetGridVelocity(CConfig *config, unsigned long iter);

	/*!
	 * \brief A virtual member.
   * \param[in] geometry - Geometry of the fine mesh.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetRestricted_GridVelocity(CGeometry *fine_mesh, CConfig *config);

	/*!
	 * \brief Find and store all vertices on a sharp corner in the geometry.
	 * \param[in] config - Definition of the particular problem.
	 */
  void ComputeSurf_Curvature(CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
  void ComputeAirfoil_Section(su2double *Plane_P0, su2double *Plane_Normal,
                              su2double MinXCoord, su2double MaxXCoord,
                              su2double MinYCoord, su2double MaxYCoord,
                              su2double MinZCoord, su2double MaxZCoord,
                              su2double *FlowVariable,
                              vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil,
                              vector<su2double> &Zcoord_Airfoil, vector<su2double> &Variable_Airfoil,
                              bool original_surface, CConfig *config);

  /*!
	 * \brief A virtual member.
	 */
  virtual su2double Compute_MaxThickness(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);
 
  /*!
	 * \brief A virtual member.
	 */
  virtual su2double Compute_Twist(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
	 * \brief A virtual member.
	 */
  virtual su2double Compute_Chord(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
	 * \brief A virtual member.
	 */
  virtual su2double Compute_Width(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
	 * \brief A virtual member.
	 */
  virtual su2double Compute_WaterLineWidth(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
	 * \brief A virtual member.
	 */
  virtual su2double Compute_Height(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_LERadius(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
	 * \brief A virtual member.
	 */
	virtual su2double Compute_Thickness(su2double *Plane_P0, su2double *Plane_Normal, su2double Location, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil, su2double &ZLoc);
	
	/*!
	 * \brief A virtual member.
	 */
	virtual su2double Compute_Area(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);
  
	/*!
	 * \brief A virtual member.
	 */
	virtual su2double Compute_Length(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
   * \brief A virtual member.
   */
  virtual void Compute_Wing_LeadingTrailing(su2double *LeadingEdge, su2double *TrailingEdge, su2double *Plane_P0, su2double *Plane_Normal, vector<su2double>
	                                       &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
   * \brief A virtual member.
   */
  virtual void Compute_Fuselage_LeadingTrailing(su2double *LeadingEdge, su2double *TrailingEdge, su2double *Plane_P0, su2double *Plane_Normal, vector<su2double>
	                                       &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_Dihedral(su2double *LeadingEdge_im1, su2double *TrailingEdge_im1,
                                     su2double *LeadingEdge_i, su2double *TrailingEdge_i);

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_Curvature(su2double *LeadingEdge_im1, su2double *TrailingEdge_im1,
                                      su2double *LeadingEdge_i, su2double *TrailingEdge_i,
                                      su2double *LeadingEdge_ip1, su2double *TrailingEdge_ip1);

  /*!
   * \brief A virtual member.
   */
  virtual void Compute_Wing(CConfig *config, bool original_surface,
                            su2double &Wing_Volume, su2double &Wing_MinMaxThickness, su2double &Wing_MaxMaxThickness, su2double &Wing_MinChord, su2double &Wing_MaxChord,
                            su2double &Wing_MinLERadius, su2double &Wing_MaxLERadius,
                            su2double &Wing_MinToC, su2double &Wing_MaxToC, su2double &Wing_ObjFun_MinToC, su2double &Wing_MaxTwist, su2double &Wing_MaxCurvature,
                            su2double &Wing_MaxDihedral);

  /*!
   * \brief A virtual member.
   */
  virtual void Compute_Fuselage(CConfig *config, bool original_surface,
  		su2double &Fuselage_Volume, su2double &Fuselage_WettedArea,
  		su2double &Fuselage_MinWidth, su2double &Fuselage_MaxWidth,
  		su2double &Fuselage_MinWaterLineWidth, su2double &Fuselage_MaxWaterLineWidth,
  		su2double &Fuselage_MinHeight, su2double &Fuselage_MaxHeight,
  		su2double &Fuselage_MaxCurvature);

  /*!
   * \brief A virtual member.
   */
  virtual void Compute_Nacelle(CConfig *config, bool original_surface,
                               su2double &Nacelle_Volume, su2double &Nacelle_MinMaxThickness, su2double &Nacelle_MaxMaxThickness,
                               su2double &Nacelle_MinChord, su2double &Nacelle_MaxChord,
                               su2double &Nacelle_MinLERadius, su2double &Nacelle_MaxLERadius,
                               su2double &Nacelle_MinToC, su2double &Nacelle_MaxToC,
                               su2double &Nacelle_ObjFun_MinToC, su2double &Nacelle_MaxTwist);

	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void FindNormal_Neighbor(CConfig *config);

  /*!
   * \brief A virtual member.
   */
  virtual void SetGlobal_to_Local_Point();

	/*!
	 * \brief A virtual member.
	 * \param[in] val_ipoint - Global point.
	 * \returns Local index that correspond with the global index.
	 */
	virtual long GetGlobal_to_Local_Point(unsigned long val_ipoint);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_ipoint - Global marker.
	 * \returns Local marker that correspond with the global index.
	 */
	virtual unsigned short GetGlobal_to_Local_Marker(unsigned short val_imarker);

  /*!
	 * \brief A virtual member.
	 * \returns Total number of nodes in a simulation across all processors (including halos).
	 */
	virtual unsigned long GetGlobal_nPoint();
  
	/*!
	 * \brief A virtual member.
	 * \returns Total number of nodes in a simulation across all processors (excluding halos).
	 */
	virtual unsigned long GetGlobal_nPointDomain();

  /*!
   * \brief A virtual member.
	 * \param[in] val_global_npoint - Global number of points in the mesh (excluding halos).
   */
  virtual void SetGlobal_nPointDomain(unsigned long val_global_npoint);

  /*!
	 * \brief A virtual member.
	 * \returns Total number of elements in a simulation across all processors.
	 */
	virtual unsigned long GetGlobal_nElem();

  /*!
   * \brief A virtual member.
   * \returns Total number of elements in a simulation across all processors (excluding halos).
   */
  virtual unsigned long GetGlobal_nElemDomain();
  
  /*!
	 * \brief A virtual member.
	 * \returns Total number of line elements in a simulation across all processors.
	 */
	virtual unsigned long GetGlobal_nElemLine();
  
  /*!
	 * \brief A virtual member.
	 * \returns Total number of triangular elements in a simulation across all processors.
	 */
	virtual unsigned long GetGlobal_nElemTria();
  
  /*!
	 * \brief A virtual member.
	 * \returns Total number of quadrilateral elements in a simulation across all processors.
	 */
	virtual unsigned long GetGlobal_nElemQuad();
  
  /*!
	 * \brief A virtual member.
	 * \returns Total number of tetrahedral elements in a simulation across all processors.
	 */
	virtual unsigned long GetGlobal_nElemTetr();
  
  /*!
	 * \brief A virtual member.
	 * \returns Total number of hexahedral elements in a simulation across all processors.
	 */
	virtual unsigned long GetGlobal_nElemHexa();
  
  /*!
	 * \brief A virtual member.
	 * \returns Total number of prism elements in a simulation across all processors.
	 */
	virtual unsigned long GetGlobal_nElemPris();
  
  /*!
	 * \brief A virtual member.
	 * \returns Total number of pyramid elements in a simulation across all processors.
	 */
	virtual unsigned long GetGlobal_nElemPyra();
  
  /*!
	 * \brief A virtual member.
	 * \return Number of line elements.
	 */
	virtual unsigned long GetnElemLine();
  
  /*!
	 * \brief A virtual member.
	 * \return Number of triangular elements.
	 */
	virtual unsigned long GetnElemTria();
  
  /*!
	 * \brief A virtual member.
	 * \return Number of quadrilateral elements.
	 */
	virtual unsigned long GetnElemQuad();
  
  /*!
	 * \brief A virtual member.
	 * \return Number of tetrahedral elements.
	 */
	virtual unsigned long GetnElemTetr();
  
  /*!
	 * \brief A virtual member.
	 * \return Number of hexahedral elements.
	 */
	virtual unsigned long GetnElemHexa();
  
  /*!
	 * \brief A virtual member.
	 * \return Number of prism elements.
	 */
	virtual unsigned long GetnElemPris();
  
  /*!
	 * \brief A virtual member.
	 * \return Number of pyramid elements.
	 */
	virtual unsigned long GetnElemPyra();

	/*!
	 * \brief Indentify geometrical planes in the mesh
	 */
	virtual void SetGeometryPlanes(CConfig *config);
  
	/*!
	 * \brief Get geometrical planes in the mesh
	 */
	virtual vector<su2double> GetGeometryPlanes();

	/*!
	 * \brief Get x coords of geometrical planes in the mesh
	 */
	virtual vector<vector<su2double> > GetXCoord();

	/*!
	 * \brief Get y coords of geometrical planes in the mesh
	 */
	virtual vector<vector<su2double> > GetYCoord();

	/*!
	 * \brief Get z coords of geometrical planes in the mesh
	 */
	virtual vector<vector<su2double> > GetZCoord();

	/*!
	 * \brief Get all points on a geometrical plane in the mesh
	 */
	virtual vector<vector<unsigned long> > GetPlanarPoints();

	/*!
	 * \brief Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f(xi), with 
	          x1 < x2 < . . . < xN , and given values yp1 and ypn for the first derivative of the interpolating
	          function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains
	          the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or
	          ypn are equal to 1 × 1030 or larger, the routine is signaled to set the corresponding boundary
	          condition for a natural spline, with zero second derivative on that boundary.
						Numerical Recipes: The Art of Scientific Computing, Third Edition in C++.
	 */
	void SetSpline(vector<su2double> &x, vector<su2double> &y, unsigned long n, su2double yp1, su2double ypn, vector<su2double> &y2);
	
	/*!
	 * \brief Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai’s in order), 
	          and given the array y2a[1..n], which is the output from spline above, and given a value of 
	          x, this routine returns a cubic-spline interpolated value y.
         	  Numerical Recipes: The Art of Scientific Computing, Third Edition in C++.
	 * \returns The interpolated value of for x.
	 */
	su2double GetSpline(vector<su2double> &xa, vector<su2double> &ya, vector<su2double> &y2a, unsigned long n, su2double x);
	  
  /*!
	 * \brief Compute the intersection between a segment and a plane.
   * \param[in] Segment_P0 - Definition of the particular problem.
	 * \param[in] Segment_P1 - Definition of the particular problem.
	 * \param[in] Plane_P0 - Definition of the particular problem.
	 * \param[in] Plane_Normal - Definition of the particular problem.
   * \param[in] Intersection - Definition of the particular problem.
   * \returns If the intersection has has been successful.
   */
  bool SegmentIntersectsPlane(su2double *Segment_P0, su2double *Segment_P1, su2double Variable_P0, su2double Variable_P1,
                              su2double *Plane_P0, su2double *Plane_Normal, su2double *Intersection, su2double &Variable_Interp);
  
  /*!
   * \brief Ray Intersects Triangle (Moller and Trumbore algorithm)
   */
  bool RayIntersectsTriangle(su2double orig[3], su2double dir[3],
                             su2double vert0[3], su2double vert1[3], su2double vert2[3],
                             su2double *intersect);
  
  /*!
   * \brief Segment Intersects Triangle
   */
  bool SegmentIntersectsTriangle(su2double point0[3], su2double point1[3],
                                 su2double vert0[3], su2double vert1[3], su2double vert2[3]);

  /*!
   * \brief Segment Intersects Line (for 2D FFD Intersection)
   */
  bool SegmentIntersectsLine(su2double point0[2], su2double point1[2], su2double vert0[2], su2double vert1[2]);

  /*!
   * \brief Register the coordinates of the mesh nodes.
   * \param[in] config
   */
  void RegisterCoordinates(CConfig *config);

  /*!
   * \brief Register the coordinates of the mesh nodes as output.
   * \param[in] config
   */
  void RegisterOutput_Coordinates(CConfig *config);

  /*!
   * \brief Update the multi-grid structure and the wall-distance.
   * \param geometry_container - Geometrical definition.
   * \param config - Config
   */
  void UpdateGeometry(CGeometry **geometry_container, CConfig *config);

  /*!
   * \brief Update the multi-grid structure for the customized boundary conditions
   * \param geometry_container - Geometrical definition.
   * \param config - Definition of the particular problem.
   */
  void UpdateCustomBoundaryConditions(CGeometry **geometry_container, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param config - Config
   */
  virtual void SetSensitivity(CConfig *config);

  /*!
   * \brief A virtual member.
   * \param config - Config
   */
  virtual void ReadUnorderedSensitivity(CConfig *config);

  /*!
   * \brief A virtual member.
   * \param iPoint - Point
   * \param iDim - Dimension
   */
  virtual su2double GetSensitivity(unsigned long iPoint, unsigned short iDim);

  /*!
   * \brief A virtual member.
   * \param iPoint - Point
   * \param iDim - Dimension
   * \param val - Value of the sensitivity
   */
  virtual void SetSensitivity(unsigned long iPoint, unsigned short iDim, su2double val);

  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker value.
	 * \param[in] val_span - span value.
	 */
  virtual su2double* GetAverageTurboNormal(unsigned short val_marker, unsigned short val_span);

  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker value.
	 * \param[in] val_span - span value.
	 */
  virtual su2double* GetAverageNormal(unsigned short val_marker, unsigned short val_span);

  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker value.
	 * \param[in] val_span - span value.
	 */
	virtual su2double GetSpanArea(unsigned short val_marker, unsigned short val_span);

  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker value.
	 * \param[in] val_span - span value.
	 */
	virtual su2double GetTurboRadius(unsigned short val_marker, unsigned short val_span);

  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker value.
	 * \param[in] val_span - span value.
	 */
	virtual su2double GetAverageTangGridVel(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   * \return The span-wise inflow tangential velocity.
	 */
	virtual su2double GetTangGridVelIn(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   * \return The span-wise outflow tangential velocity.
	 */
	virtual su2double GetTangGridVelOut(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   * \return The span-wise inflow area.
	 */
	virtual su2double GetSpanAreaIn(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
	 * \return The span-wise outflow area.
	 */
	virtual su2double GetSpanAreaOut(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   * \return The span-wise inflow radius.
	 */
	virtual su2double GetTurboRadiusIn(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   * \return The span-wise outflow radius.
	 */
	virtual su2double GetTurboRadiusOut(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   */
	virtual void SetTangGridVelIn(su2double value, unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   */
	virtual void SetTangGridVelOut(su2double value, unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   */
	virtual void SetSpanAreaIn(su2double value, unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   */
	virtual void SetSpanAreaOut(su2double value, unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   */
	virtual void SetTurboRadiusIn(su2double value, unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   */
	virtual void SetTurboRadiusOut(su2double value, unsigned short val_marker, unsigned short val_span);

  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - marker value.
	 * \param[in] val_span - span value.
	 */
	virtual unsigned long GetnTotVertexSpan(unsigned short val_marker, unsigned short val_span);

/*!
 * \brief A virtual member.
 * \param[in] val_marker - marker value.
 * \param[in] val_span - span value.
 */
  virtual su2double GetMinAngularCoord(unsigned short val_marker, unsigned short val_span);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - marker value.
   * \param[in] val_span - span value.
   */
  virtual su2double GetMaxAngularCoord(unsigned short val_marker, unsigned short val_span);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - marker value.
   * \param[in] val_span - span value.
   */
  virtual su2double GetMinRelAngularCoord(unsigned short val_marker, unsigned short val_span);
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - marker value.
   * \param[in] val_span - span value.
   */
  virtual su2double* GetAverageGridVel(unsigned short val_marker, unsigned short val_span);

  /*!
   * \brief A virtual member.
   * \param config - Config
   */
  virtual void Check_Periodicity(CConfig *config);

  /*!
   * \brief Get the value of the customized temperature at a specified vertex on a specified marker.
   * \param[in] val_marker - Marker value
   * \param[in] val_vertex - Boundary vertex value
   */
  su2double GetCustomBoundaryTemperature(unsigned short val_marker, unsigned long val_vertex);

  /*!
   * \brief Set the value of the customized temperature at a specified vertex on a specified marker.
   * \param[in] val_marker - Marker value
   * \param[in] val_vertex - Boundary vertex value
   * \param[in] val_customBoundaryTemperature - Value of the temperature.
   */
  void SetCustomBoundaryTemperature(unsigned short val_marker, unsigned long val_vertex, su2double val_customBoundaryTemperature);

  /*!
   * \brief Get the value of the customized normal heat flux at a specified vertex on a specified marker.
   * \param[in] val_marker - Marker value
   * \param[in] val_vertex - Boundary vertex value
   */
  su2double GetCustomBoundaryHeatFlux(unsigned short val_marker, unsigned long val_vertex);

  /*!
   * \brief Set the value of the customized normal heat flux at a specified vertex on a specified marker.
   * \param[in] val_marker - Marker value
   * \param[in] val_vertex - Boundary vertex value
   * \param[in] val_customBoundaryHeatFlux - Value of the normal heat flux.
   */
  void SetCustomBoundaryHeatFlux(unsigned short val_marker, unsigned long val_vertex, su2double val_customBoundaryHeatFlux);
  
  /*!
   * \brief Filter values given at the element CG by performing a weighted average over a radial neighbourhood.
   * \param[in] filter_radius - Parameter defining the size of the neighbourhood.
   * \param[in] kernels - Kernel types and respective parameter, size of vector defines number of filter recursions.
   * \param[in] input_values - "Raw" values.
   * \param[out] output_values - Filtered values.
   */
  void FilterValuesAtElementCG(const vector<su2double> filter_radius, const vector<pair<unsigned short,su2double> > &kernels,
                               const su2double *input_values, su2double *output_values) const;
  
  /*!
   * \brief Build the global (entire mesh!) adjacency matrix for the elements in compressed format.
   *        Used by FilterValuesAtElementCG to search for geometrically close neighbours.
   * \param[out] neighbour_start - i'th position stores the start position in "neighbour_idx" for the immediate
   *             neighbours of global element "i". Size nElemDomain+1
   * \param[out] neighbour_idx - Global index of the neighbours, mush be NULL on entry and free'd by calling function.
   */
  void GetGlobalElementAdjacencyMatrix(vector<unsigned long> &neighbour_start, long *&neighbour_idx) const;
  
  /*!
   * \brief Get the neighbours of the global element in the first position of "neighbours" that are within "radius" of it.
   * \param[in] iElem_global - Element of interest.
   * \param[in] radius - Parameter defining the size of the neighbourhood.
   * \param[in] neighbour_start - See GetGlobalElementAdjacencyMatrix.
   * \param[in] neighbour_idx - See GetGlobalElementAdjacencyMatrix.
   * \param[in] cg_elem - Global element centroid coordinates in row major format {x0,y0,x1,y1,...}. Size nDim*nElemDomain.
   * \param[in,out] neighbours - The neighbours of iElem_global.
   */
  void GetRadialNeighbourhood(const unsigned long iElem_global, const passivedouble radius,
                              const vector<unsigned long> &neighbour_start, const long *neighbour_idx,
                              const su2double *cg_elem, vector<long> &neighbours) const;

  /*!
   * \brief Compute and store the volume of the elements.
   * \param[in] config - Problem configuration.
   */
  void SetElemVolume(CConfig *config);
  
  /*!
   * \brief Set the multigrid index for the current geometry object.
   * \param[in] val_iMesh - Multigrid index for current geometry object.
   */
  void SetMGLevel(unsigned short val_iMesh);
  
  /*!
   * \brief Get the multigrid index for the current geometry object.
   * \return Multigrid index for current geometry object.
   */
  unsigned short GetMGLevel(void);

};

/*!
 * \class CPhysicalGeometry
 * \brief Class for reading a defining the primal grid which is read from the
 *        grid file in .su2 or .cgns format.
 * \author F. Palacios, T. Economon, J. Alonso
 */
class CPhysicalGeometry : public CGeometry {

  map<unsigned long, unsigned long> Global_to_Local_Point; /*!< \brief Global-local indexation for the points. */
  long *Local_to_Global_Point;				/*!< \brief Local-global indexation for the points. */
  unsigned short *Local_to_Global_Marker;	/*!< \brief Local to Global marker. */
  unsigned short *Global_to_Local_Marker;	/*!< \brief Global to Local marker. */
  unsigned long *adj_counter; /*!< \brief Adjacency counter. */
  unsigned long **adjacent_elem; /*!< \brief Adjacency element list. */
  su2double* Sensitivity; /*! <\brief Vector holding the sensitivities at each point. */

  vector<vector<unsigned long> > Neighbors;
  map<unsigned long, unsigned long> Color_List;
  vector<string> Marker_Tags;
  unsigned long nLocal_Point,
  nLocal_PointDomain,
  nLocal_PointGhost,
  nLocal_PointPeriodic,
  nLocal_Elem,
  nLocal_Bound_Elem,
  nGlobal_Elem,
  nGlobal_Bound_Elem,
  nLocal_Line,
  nLocal_BoundTria,
  nLocal_BoundQuad,
  nLinear_Line,
  nLinear_BoundTria,
  nLinear_BoundQuad,
  nLocal_Tria,
  nLocal_Quad,
  nLocal_Tetr,
  nLocal_Hexa,
  nLocal_Pris,
  nLocal_Pyra;
  unsigned long nMarker_Global;
  su2double *Local_Coords;
  unsigned long *Local_Points;
  unsigned long *Local_Colors;
  unsigned long *Conn_Line;
  unsigned long *Conn_BoundTria;
  unsigned long *Conn_BoundQuad;
  unsigned long *Conn_Line_Linear;
  unsigned long *Conn_BoundTria_Linear;
  unsigned long *Conn_BoundQuad_Linear;
  unsigned long *Conn_Tria;
  unsigned long *Conn_Quad;
  unsigned long *Conn_Tetr;
  unsigned long *Conn_Hexa;
  unsigned long *Conn_Pris;
  unsigned long *Conn_Pyra;
  unsigned long *ID_Line;
  unsigned long *ID_BoundTria;
  unsigned long *ID_BoundQuad;
  unsigned long *ID_Line_Linear;
  unsigned long *ID_BoundTria_Linear;
  unsigned long *ID_BoundQuad_Linear;
  unsigned long *ID_Tria;
  unsigned long *ID_Quad;
  unsigned long *ID_Tetr;
  unsigned long *ID_Hexa;
  unsigned long *ID_Pris;
  unsigned long *ID_Pyra;
  unsigned long *Elem_ID_Line;
  unsigned long *Elem_ID_BoundTria;
  unsigned long *Elem_ID_BoundQuad;
  unsigned long *Elem_ID_Line_Linear;
  unsigned long *Elem_ID_BoundTria_Linear;
  unsigned long *Elem_ID_BoundQuad_Linear;

public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CPhysicalGeometry(void);

	/*! 
	 * \overload
	 * \brief Reads the geometry of the grid and adjust the boundary 
	 *        conditions with the configuration file. 
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_mesh_filename - Name of the file with the grid information.
	 * \param[in] val_format - Format of the file with the grid information.
	 * \param[in] val_iZone - Domain to be read from the grid file.
	 * \param[in] val_nZone - Total number of domains in the grid file.
	 */
	CPhysicalGeometry(CConfig *config, unsigned short val_iZone, unsigned short val_nZone);
  
  /*!
	 * \overload
	 * \brief Accepts a geometry container holding a linearly partitioned grid
   *        with coloring performed by ParMETIS, and this routine distributes
   *        the points and cells to all partitions based on the coloring.
   * \param[in] geometry - Definition of the geometry container holding the initial linear partitions of the grid + coloring.
	 * \param[in] config - Definition of the particular problem.
	 */
  CPhysicalGeometry(CGeometry *geometry, CConfig *config);

  /*!
   * \overload
   * \brief Accepts a geometry container holding a linearly partitioned grid
   *        with coloring performed by ParMETIS, and this routine distributes
   *        the points and cells to all partitions based on the coloring.
   * \param[in] geometry - Definition of the geometry container holding the initial linear partitions of the grid + coloring.
   * \param[in] config - Definition of the particular problem.
   */
  CPhysicalGeometry(CGeometry *geometry, CConfig *config, bool val_flag);

	/*!
	 * \brief Destructor of the class.
	 */
	~CPhysicalGeometry(void);

  /*!
   * \brief Distributes the coloring from ParMETIS so that each rank has complete information about the local grid points.
   * \param[in] geometry - Definition of the geometry container holding the initial linear partitions of the grid + coloring.
   * \param[in] config - Definition of the particular problem.
   */
  void DistributeColoring(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Distribute the grid points, including ghost points, across all ranks based on a ParMETIS coloring.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void DistributePoints(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Distribute the connectivity for a single volume element type across all ranks based on a ParMETIS coloring.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being distributed.
   */
  void DistributeVolumeConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);

  /*!
   * \brief Distribute the connectivity for a single surface element type in all markers across all ranks based on a ParMETIS coloring.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being distributed.
   */
  void DistributeSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);

  /*!
   * \brief Broadcast the marker tags for all boundaries from the master rank to all other ranks.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void DistributeMarkerTags(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Partition the marker connectivity held on the master rank according to a linear partitioning.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being distributed.
   */
  void PartitionSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);

  /*!
   * \brief Load the local grid points after partitioning (owned and ghost) into the geometry class objects.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void LoadPoints(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Load the local volume elements after partitioning (owned and ghost) into the geometry class objects.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void LoadVolumeElements(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Load the local surface elements after partitioning (owned and ghost) into the geometry class objects.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void LoadSurfaceElements(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Routine to launch non-blocking sends and recvs amongst all processors.
   * \param[in] bufSend - Buffer of data to be sent.
   * \param[in] nElemSend - Array containing the number of elements to send to other processors in cumulative storage format.
   * \param[in] sendReq - Array of MPI send requests.
   * \param[in] bufRecv - Buffer of data to be received.
   * \param[in] nElemSend - Array containing the number of elements to receive from other processors in cumulative storage format.
   * \param[in] sendReq - Array of MPI recv requests.
   * \param[in] countPerElem - Pieces of data per element communicated.
   */
  void InitiateCommsAll(void *bufSend,
                        int *nElemSend,
                        SU2_MPI::Request *sendReq,
                        void *bufRecv,
                        int *nElemRecv,
                        SU2_MPI::Request *recvReq,
                        unsigned short countPerElem,
                        unsigned short commType);

  /*!
   * \brief Routine to complete the set of non-blocking communications launched with InitiateComms() with MPI_Waitany().
   * \param[in] nSends - Number of sends to be completed.
   * \param[in] sendReq - Array of MPI send requests.
   * \param[in] nRecvs - Number of receives to be completed.
   * \param[in] sendReq - Array of MPI recv requests.
   */
  void CompleteCommsAll(int nSends,
                        SU2_MPI::Request *sendReq,
                        int nRecvs,
                        SU2_MPI::Request *recvReq);

  /*!
	 * \brief Set the send receive boundaries of the grid.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_domain - Number of domains for parallelization purposes.
	 */
	void SetSendReceive(CConfig *config);
  
  /*!
	 * \brief Set the send receive boundaries of the grid.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_domain - Number of domains for parallelization purposes.
	 */
	void SetBoundaries(CConfig *config);

  /*!
   * \brief Set the local index that correspond with the global numbering index.
   */
  void SetGlobal_to_Local_Point();

	/*!
	 * \brief Get the local index that correspond with the global numbering index.
	 * \param[in] val_ipoint - Global point.
	 * \returns Local index that correspond with the global index, -1 if not found on the current rank (process).
	 */
	long GetGlobal_to_Local_Point(unsigned long val_ipoint);
  
	/*!
	 * \brief Get the local marker that correspond with the global marker.
	 * \param[in] val_ipoint - Global marker.
	 * \returns Local marker that correspond with the global index.
	 */
	unsigned short GetGlobal_to_Local_Marker(unsigned short val_imarker);
  
  /*!
   * \brief Reads the geometry of the grid and adjust the boundary
   *        conditions with the configuration file in parallel (for parmetis).
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_mesh_filename - Name of the file with the grid information.
   * \param[in] val_format - Format of the file with the grid information.
   * \param[in] val_iZone - Domain to be read from the grid file.
   * \param[in] val_nZone - Total number of domains in the grid file.
   */
  void Read_SU2_Format_Parallel(CConfig *config, string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone);

  /*!
   * \brief Reads the geometry of the grid and adjust the boundary
   *        conditions with the configuration file in parallel (for parmetis).
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_mesh_filename - Name of the file with the grid information.
   * \param[in] val_format - Format of the file with the grid information.
   * \param[in] val_iZone - Domain to be read from the grid file.
   * \param[in] val_nZone - Total number of domains in the grid file.
   */
  void Read_CGNS_Format_Parallel(CConfig *config, string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone);

  /*!
   * \brief Reads for the FEM solver the geometry of the grid and adjust the boundary
   *        conditions with the configuration file in parallel (for parmetis).
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_mesh_filename - Name of the file with the grid information.
   * \param[in] val_iZone - Domain to be read from the grid file.
   * \param[in] val_nZone - Total number of domains in the grid file.
   */
  void Read_SU2_Format_Parallel_FEM(CConfig *config, string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone);

  /*!
   * \brief Reads for the FEM solver the geometry of the grid and adjust the boundary
   *        conditions with the configuration file in parallel (for parmetis).
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_mesh_filename - Name of the file with the grid information.
   * \param[in] val_iZone - Domain to be read from the grid file.
   * \param[in] val_nZone - Total number of domains in the grid file.
   */
  void Read_CGNS_Format_Parallel_FEM(CConfig *config, string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone);

	/*!
	 * \brief Find repeated nodes between two elements to identify the common face.
	 * \param[in] first_elem - Identification of the first element.
	 * \param[in] second_elem - Identification of the second element.
	 * \param[in] face_first_elem - Index of the common face for the first element.
	 * \param[in] face_second_elem - Index of the common face for the second element.		 
	 * \return It provides 0 or 1 depending if there is a common face or not.
	 */
	bool FindFace(unsigned long first_elem, unsigned long second_elem, unsigned short &face_first_elem, 
			unsigned short &face_second_elem);

	/*! 
	 * \brief Computes the distance to the nearest no-slip wall for each grid node.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeWall_Distance(CConfig *config);

	/*! 
	 * \brief Compute surface area (positive z-direction) for force coefficient non-dimensionalization.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPositive_ZArea(CConfig *config);
  
	/*! 
	 * \brief Set points which surround a point.
	 */
	void SetPoint_Connectivity(void);
  
  /*!
	 * \brief Set a renumbering using a Reverse Cuthill-McKee Algorithm
   * \param[in] config - Definition of the particular problem.
	 */
	void SetRCM_Ordering(CConfig *config);
  
	/*!
	 * \brief Function declaration to avoid partially overridden classes.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */
	void SetPoint_Connectivity(CGeometry *geometry);

	/*! 
	 * \brief Set elements which surround an element.
	 */
	void SetElement_Connectivity(void);

	/*! 
	 * \brief Set the volume element associated to each boundary element.
	 */
	void SetBoundVolume(void);

	/*! 
	 * \brief Set boundary vertex.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetVertex(CConfig *config);

	/*!
	 * \brief Set number of span wise level for turbomachinery computation.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeNSpan(CConfig *config, unsigned short val_iZone, unsigned short marker_flag, bool allocate);

	/*!
	 * \brief Set turbo boundary vertex.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetTurboVertex(CConfig *config,unsigned short val_iZone, unsigned short marker_flag, bool allocate);

/*!
 * \brief update turbo boundary vertex.
 * \param[in] config - Definition of the particular problem.
 */
void UpdateTurboVertex(CConfig *config,unsigned short val_iZone, unsigned short marker_flag);

	/*!
	 * \brief Set turbo boundary vertex.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetAvgTurboValue(CConfig *config, unsigned short val_iZone, unsigned short marker_flag, bool allocate);

	/*!
	 * \brief Set turbo boundary vertex.
	 * \param[in] config - Definition of the particular problem.
	 */
	void GatherInOutAverageValues(CConfig *config, bool allocate);

	/*! 
	 * \brief Set the center of gravity of the face, elements and edges.
	 */
	void SetCoord_CG(void);

	/*! 
	 * \brief Set the edge structure of the control volume.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] action - Allocate or not the new elements.
	 */
	void SetControlVolume(CConfig *config, unsigned short action);

	/*!
	 * \brief Visualize the structure of the control volume(s).
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] action - Allocate or not the new elements.
	 */
  void VisualizeControlVolume(CConfig *config, unsigned short action);
  
	/*!
	 * \brief Mach the near field boundary condition.
	 * \param[in] config - Definition of the particular problem.
	 */
	void MatchNearField(CConfig *config);
  
  /*!
	 * \brief Mach the near field boundary condition.
	 * \param[in] config - Definition of the particular problem.
	 */
	void MatchActuator_Disk(CConfig *config);

	/*! 
	 * \brief Mach the interface boundary condition.
	 * \param[in] config - Definition of the particular problem.
	 */
	void MatchInterface(CConfig *config);
  
  /*!
   * \brief Mach the periodic boundary conditions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_periodic - Index of the first periodic face in a pair.
   */
  void MatchPeriodic(CConfig *config, unsigned short val_periodic);

	/*! 
	 * \brief Set boundary vertex structure of the control volume.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] action - Allocate or not the new elements.
	 */
	void SetBoundControlVolume(CConfig *config, unsigned short action);

  /*!
   * \brief Set the maximum cell-center to cell-center distance for CVs.
   * \param[in] config - Definition of the particular problem.
   */
  void SetMaxLength(CConfig* config);

	/*! 
	 * \brief Set the Tecplot file.
	 * \param[in] config_filename - Name of the file where the Tecplot 
	 *            information is going to be stored.
   * \param[in] new_file - Create a new file.
	 */
	void SetTecPlot(char config_filename[MAX_STRING_SIZE], bool new_file);

	/*! 
	 * \brief Set the output file for boundaries in Tecplot
	 * \param[in] config - Definition of the particular problem.		 
	 * \param[in] mesh_filename - Name of the file where the Tecplot 
	 *            information is going to be stored.   
   * \param[in] new_file - Create a new file.
	 */
	void SetBoundTecPlot(char mesh_filename[MAX_STRING_SIZE], bool new_file, CConfig *config);

	/*! 
	 * \brief Check the volume element orientation.
	 * \param[in] config - Definition of the particular problem.		 
	 */
	void Check_IntElem_Orientation(CConfig *config);
  
  /*!
	 * \brief Check the volume element orientation.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Check_BoundElem_Orientation(CConfig *config);

	/*! 
	 * \brief Set the domains for grid grid partitioning using METIS.
	 * \param[in] config - Definition of the particular problem.		 
	 */
	void SetColorGrid(CConfig *config);
  
  /*!
   * \brief Set the domains for grid grid partitioning using ParMETIS.
   * \param[in] config - Definition of the particular problem.
   */
  void SetColorGrid_Parallel(CConfig *config);

  /*!
   * \brief Set the domains for FEM grid partitioning using ParMETIS.
   * \param[in] config - Definition of the particular problem.
   */
  void SetColorFEMGrid_Parallel(CConfig *config);

  /*!
   * \brief Compute the weights of the FEM graph for ParMETIS.
   * \param[in]  config                       - Definition of the particular problem.
   * \param[in]  localFaces                   - Vector, which contains the element faces of this rank.
   * \param[in]  adjacency                    - Neighbors of the element.
   * \param[in]  mapExternalElemIDToTimeLevel - Map from the external element ID's to their time level
                                                and number of DOFs.
   * \param[out] vwgt                         - Weights of the vertices of the graph, i.e. the elements.
   * \param[out] adjwgt                       - Weights of the edges of the graph.
   */
  void ComputeFEMGraphWeights(
              CConfig                                    *config,
              const vector<CFaceOfElement>               &localFaces,
              const vector<vector<unsigned long> >       &adjacency,
              const map<unsigned long, CUnsignedShort2T> &mapExternalElemIDToTimeLevel,
                    vector<su2double>                    &vwgt,
                    vector<vector<su2double> >           &adjwgt);

  /*!
   * \brief Determine the donor elements for the boundary elements on viscous
            wall boundaries when wall functions are used.
   * \param[in]  config - Definition of the particular problem.
   */
  void DetermineDonorElementsWallFunctions(CConfig *config);

  /*!
   * \brief Determine whether or not the Jacobians of the elements and faces
            are constant and a length scale of the elements.
   * \param[in]  config - Definition of the particular problem.
   */
  void DetermineFEMConstantJacobiansAndLenScale(CConfig *config);

  /*!
   * \brief Determine the neighboring information for periodic faces of a FEM grid.
   * \param[in]     config      - Definition of the particular problem.
   * \param[in,out] localFaces  - Vector, which contains the element faces of this rank.
   */
  void DeterminePeriodicFacesFEMGrid(CConfig                *config,
                                     vector<CFaceOfElement> &localFaces);

  /*!
   * \brief Determine the time level of the elements when time accurate
            local time stepping is employed.
   * \param[in]  config                       - Definition of the particular problem.
   * \param[in]  localFaces                   - Vector, which contains the element
                                                faces of this rank.
   * \param[out] mapExternalElemIDToTimeLevel - Map from the external element ID's to
                                                their time level and number of DOFs.
   */
  void DetermineTimeLevelElements(CConfig                              *config,
                                  const vector<CFaceOfElement>         &localFaces,
                                  map<unsigned long, CUnsignedShort2T> &mapExternalElemIDToTimeLevel);

  /*!
   * \brief Set the rotational velocity at each node.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iZone - Index of the current zone.
   */
  void SetRotationalVelocity(CConfig *config, unsigned short val_iZone, bool print);

  /*!
   * \brief Set the rotational velocity of the points on the shroud markers to 0.
   * \param[in] config - Definition of the particular problem.
   */
  void SetShroudVelocity(CConfig *config);

  /*!
   * \brief Set the translational velocity at each node.
   * \param[in] config - Definition of the particular problem.
   */
  void SetTranslationalVelocity(CConfig *config, unsigned short val_iZone, bool print);

  /*!
   * \brief Set the grid velocity via finite differencing at each node.
   * \param[in] config - Definition of the particular problem.
   */
  void SetGridVelocity(CConfig *config, unsigned long iter);

	/*! 
	 * \brief Do an implicit smoothing of the grid coordinates.
	 * \param[in] val_nSmooth - Number of smoothing iterations.
	 * \param[in] val_smooth_coeff - Relaxation factor.
	 * \param[in] config - Definition of the particular problem.		 
	 */	
	void SetCoord_Smoothing(unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config);

	/*! 
	 * \brief Write the .su2 file.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_mesh_out_filename - Name of the output file.
	 */	
	void SetMeshFile(CConfig *config, string val_mesh_out_filename);

	/*! 
	 * \brief Compute some parameters about the grid quality.
	 * \param[out] statistics - Information about the grid quality, statistics[0] = (r/R)_min, statistics[1] = (r/R)_ave.		 
	 */	
	void GetQualityStatistics(su2double *statistics);

	/*!
	 * \brief Find and store all vertices on a sharp corner in the geometry.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeSurf_Curvature(CConfig *config);

	/*! 
	 * \brief Find and store the closest neighbor to a vertex.
	 * \param[in] config - Definition of the particular problem.
	 */
	void FindNormal_Neighbor(CConfig *config);

  /*!
	 * \brief Retrieve total number of nodes in a simulation across all processors (including halos).
	 * \returns Total number of nodes in a simulation across all processors (including halos).
	 */
	unsigned long GetGlobal_nPoint();
  
	/*!
	 * \brief Retrieve total number of nodes in a simulation across all processors (excluding halos).
	 * \returns Total number of nodes in a simulation across all processors (excluding halos).
	 */
	unsigned long GetGlobal_nPointDomain();
  
  /*!
	 * \brief Retrieve total number of elements in a simulation across all processors.
	 * \returns Total number of elements in a simulation across all processors.
	 */
  unsigned long GetGlobal_nElem();
  
  /*!
   * \brief  Retrieve total number of elements in a simulation across all processors (excluding halos).
   * \returns Total number of elements in a simulation across all processors (excluding halos).
   */
  unsigned long GetGlobal_nElemDomain();

  /*!
	 * \brief Retrieve total number of triangular elements in a simulation across all processors.
	 * \returns Total number of line elements in a simulation across all processors.
	 */
	unsigned long GetGlobal_nElemLine();
  
  /*!
	 * \brief Retrieve total number of triangular elements in a simulation across all processors.
	 * \returns Total number of triangular elements in a simulation across all processors.
	 */
	unsigned long GetGlobal_nElemTria();
  
  /*!
	 * \brief Retrieve total number of quadrilateral elements in a simulation across all processors.
	 * \returns Total number of quadrilateral elements in a simulation across all processors.
	 */
	unsigned long GetGlobal_nElemQuad();
  
  /*!
	 * \brief Retrieve total number of tetrahedral elements in a simulation across all processors.
	 * \returns Total number of tetrahedral elements in a simulation across all processors.
	 */
	unsigned long GetGlobal_nElemTetr();
  
  /*!
	 * \brief Retrieve total number of hexahedral elements in a simulation across all processors.
	 * \returns Total number of hexahedral elements in a simulation across all processors.
	 */
	unsigned long GetGlobal_nElemHexa();
  
  /*!
	 * \brief Retrieve total number of prism elements in a simulation across all processors.
	 * \returns Total number of prism elements in a simulation across all processors.
	 */
	unsigned long GetGlobal_nElemPris();
  
  /*!
	 * \brief Retrieve total number of pyramid elements in a simulation across all processors.
	 * \returns Total number of pyramid elements in a simulation across all processors.
	 */
	unsigned long GetGlobal_nElemPyra();
  
  /*!
	 * \brief Get number of triangular elements.
	 * \return Number of line elements.
	 */
	unsigned long GetnElemLine();
  
  /*!
	 * \brief Get number of triangular elements.
	 * \return Number of triangular elements.
	 */
	unsigned long GetnElemTria();
  
  /*!
	 * \brief Get number of quadrilateral elements.
	 * \return Number of quadrilateral elements.
	 */
	unsigned long GetnElemQuad();
  
  /*!
	 * \brief Get number of tetrahedral elements.
	 * \return Number of tetrahedral elements.
	 */
	unsigned long GetnElemTetr();
  
  /*!
	 * \brief Get number of hexahedral elements.
	 * \return Number of hexahedral elements.
	 */
	unsigned long GetnElemHexa();
  
  /*!
	 * \brief Get number of prism elements.
	 * \return Number of prism elements.
	 */
	unsigned long GetnElemPris();
  
  /*!
	 * \brief Get number of pyramid elements.
	 * \return Number of pyramid elements.
	 */
	unsigned long GetnElemPyra();

	/*!
	 * \brief Indentify geometrical planes in the mesh
	 */
	void SetGeometryPlanes(CConfig *config);
  
	/*!
	 * \brief Get geometrical planes in the mesh
	 */
	vector<su2double> GetGeometryPlanes();

	/*!
	 * \brief Get x coords of geometrical planes in the mesh
	 */
	vector<vector<su2double> > GetXCoord();

	/*!
	 * \brief Get y coords of geometrical planes in the mesh
	 */
	vector<vector<su2double> > GetYCoord();

	/*!
	 * \brief Get z coords of geometrical planes in the mesh
	 */
	vector<vector<su2double> > GetZCoord();

	/*!
	 * \brief Get all points on a geometrical plane in the mesh
	 */
	vector<vector<unsigned long> > GetPlanarPoints();
  
  /*!
   * \brief Read the sensitivity from an input file.
   * \param[in] config - Definition of the particular problem.
   */
  void SetBoundSensitivity(CConfig *config);
  
  /*!
   * \brief Compute the maximum thickness of an airfoil.
   * \returns Maximum thickness at a particular seccion.
   */
  su2double Compute_MaxThickness(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);
  
  /*!
   * \brief Compute the twist of an airfoil.
   * \returns Twist at a particular seccion.
   */
  su2double Compute_Twist(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
   * \brief Compute the leading/trailing edge location of an airfoil.
   */
  void Compute_Wing_LeadingTrailing(su2double *LeadingEdge, su2double *TrailingEdge, su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil,
                               vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
    * \brief Compute the leading/trailing edge location of a fuselage.
    */
   void Compute_Fuselage_LeadingTrailing(su2double *LeadingEdge, su2double *TrailingEdge, su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil,
                                vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
   * \brief Compute the chord of an airfoil.
   * \returns Chord of an airfoil.
   */
  su2double Compute_Chord(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
   * \brief Compute the chord of an airfoil.
   * \returns Chord of an airfoil.
   */
  su2double Compute_Width(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
   * \brief Compute the chord of an airfoil.
   * \returns Chord of an airfoil.
   */
  su2double Compute_WaterLineWidth(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
   * \brief Compute the chord of an airfoil.
   * \returns Chord of an airfoil.
   */
  su2double Compute_Height(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
   * \brief Compute the chord of an airfoil.
   * \returns Chord of an airfoil.
   */
  su2double Compute_LERadius(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
   * \brief Compute the thickness of an airfoil.
   */
  su2double Compute_Thickness(su2double *Plane_P0, su2double *Plane_Normal, su2double Location, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil, su2double &ZLoc);

  /*!
   * \brief Compute the area of an airfoil.
   * \returns Area of an airfoil.
   */
  su2double Compute_Area(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
   * \brief Compute the length of an airfoil.
   * \returns Area of an airfoil.
   */
  su2double Compute_Length(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil);

  /*!
   * \brief Compute the dihedral of a wing.
   * \returns Dihedral at a particular seccion.
   */
  su2double Compute_Dihedral(su2double *LeadingEdge_im1, su2double *TrailingEdge_im1,
                             su2double *LeadingEdge_i, su2double *TrailingEdge_i);

  /*!
   * \brief Compute the curvature of a wing.
   */
  su2double Compute_Curvature(su2double *LeadingEdge_im1, su2double *TrailingEdge_im1,
                              su2double *LeadingEdge_i, su2double *TrailingEdge_i,
                              su2double *LeadingEdge_ip1, su2double *TrailingEdge_ip1);
  
  /*!
   * \brief Evaluate geometrical parameters of a wing.
   */
  void Compute_Wing(CConfig *config, bool original_surface,
                    su2double &Wing_Volume, su2double &Wing_MinMaxThickness, su2double &Wing_MaxMaxThickness,
                    su2double &Wing_MinChord, su2double &Wing_MaxChord,
                    su2double &Wing_MinLERadius, su2double &Wing_MaxLERadius,
                    su2double &Wing_MinToC, su2double &Wing_MaxToC,
                    su2double &Wing_ObjFun_MinToC, su2double &Wing_MaxTwist,
                    su2double &Wing_MaxCurvature, su2double &Wing_MaxDihedral);

  /*!
   * \brief Evaluate geometrical parameters of a wing.
   */
  void Compute_Fuselage(CConfig *config, bool original_surface,
  		su2double &Fuselage_Volume, su2double &Fuselage_WettedArea,
  		su2double &Fuselage_MinWidth, su2double &Fuselage_MaxWidth,
  		su2double &Fuselage_MinWaterLineWidth, su2double &Fuselage_MaxWaterLineWidth,
  		su2double &Fuselage_MinHeight, su2double &Fuselage_MaxHeight,
  		su2double &Fuselage_MaxCurvature);
  
  /*!
   * \brief Evaluate geometrical parameters of a wing.
   */
  void Compute_Nacelle(CConfig *config, bool original_surface,
                       su2double &Nacelle_Volume, su2double &Nacelle_MinMaxThickness, su2double &Nacelle_MaxMaxThickness,
                       su2double &Nacelle_MinChord, su2double &Nacelle_MaxChord,
                       su2double &Nacelle_MinLERadius, su2double &Nacelle_MaxLERadius,
                       su2double &Nacelle_MinToC, su2double &Nacelle_MaxToC,
                       su2double &Nacelle_ObjFun_MinToC, su2double &Nacelle_MaxTwist);

  /*!
   * \brief Read the sensitivity from adjoint solution file and store it.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSensitivity(CConfig *config);

  /*!
   * \brief Read the sensitivity from unordered ASCII adjoint solution file and store it.
   * \param[in] config - Definition of the particular problem.
   */
  void ReadUnorderedSensitivity(CConfig *config);

  /*!
   * \brief Get the Sensitivity at a specific point.
   * \param[in] iPoint - The point where to get the sensitivity.
   * \param[in] iDim - The component of the dim. vector.
   * \returns The sensitivity at point iPoint and dim. iDim.
   */
  su2double GetSensitivity(unsigned long iPoint, unsigned short iDim);

  /*!
   * \brief Set the Sensitivity at a specific point.
   * \param[in] iPoint - The point where to get the sensitivity.
   * \param[in] iDim - The component of the dim. vector.
   * \param[in] val - Value of the sensitivity.
   */
  void SetSensitivity(unsigned long iPoint, unsigned short iDim, su2double val);
  
  /*!
   * \brief Check the mesh for periodicity and deactivate multigrid if periodicity is found.
   * \param[in] config - Definition of the particular problem.
   */
  void Check_Periodicity(CConfig *config);

  /*!
	 * \brief Get the average normal at a specific span for a given marker in the turbomachinery reference of frame.
	 * \param[in] val_marker - marker value.
	 * \param[in] val_span - span value.
   * \return The span-wise averaged turbo normal.
	 */
  su2double* GetAverageTurboNormal(unsigned short val_marker, unsigned short val_span);

  /*!
	 * \brief Get the average normal at a specific span for a given marker.
	 * \param[in] val_marker - marker value.
	 * \param[in] val_span - span value.
   * \return The span-wise averaged normal.
	 */
  su2double* GetAverageNormal(unsigned short val_marker, unsigned short val_span);

  /*!
	 * \brief Get the value of the total area for each span.
	 * \param[in] val_marker - marker value.
	 * \param[in] val_span - span value.
   * \return The span-wise area.
	 */
	su2double GetSpanArea(unsigned short val_marker, unsigned short val_span);

  /*!
	 * \brief Get the value of the total area for each span.
	 * \param[in] val_marker - marker value.
	 * \param[in] val_span - span value.
   * \return The span-wise averaged turbo normal.
	 */
	su2double GetTurboRadius(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Get the value of the average tangential rotational velocity for each span.
	 * \param[in] val_marker - marker value.
	 * \param[in] val_span - span value.
   * \return The span-wise averaged tangential velocity.
	 */
	su2double GetAverageTangGridVel(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Get the value of the inflow tangential velocity at each span.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   * \return The span-wise inflow tangential velocity.
	 */
	su2double GetTangGridVelIn(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Get the value of the outflow tangential velocity at each span.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   * \return The span-wise outflow tangential velocity.
	 */
	su2double GetTangGridVelOut(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Get the value of the inflow area at each span.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   * \return The span-wise inflow area.
	 */
	su2double GetSpanAreaIn(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Get the value of the outflow area at each span.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
	 * \return The span-wise outflow area.
	 */
	su2double GetSpanAreaOut(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Get the value of the inflow radius at each span.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   * \return The span-wise inflow radius.
	 */
	su2double GetTurboRadiusIn(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Get the value of the outflow radius at each span.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   * \return The span-wise outflow radius.
	 */
	su2double GetTurboRadiusOut(unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Set the value of the inflow tangential velocity at each span.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   */
	void SetTangGridVelIn(su2double value, unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Set the value of the outflow tangential velocity at each span.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   */
	void SetTangGridVelOut(su2double value, unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Get the value of the inflow area at each span.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   */
	void SetSpanAreaIn(su2double value, unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Set the value of the outflow area at each span.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   */
	void SetSpanAreaOut(su2double value, unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Set the value of the inflow radius at each span.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   */
	void SetTurboRadiusIn(su2double value, unsigned short val_marker, unsigned short val_span);

	/*!
	 * \brief Set the value of the outflow radius at each span.
	 * \param[in] val_marker - marker turbo-performance value.
	 * \param[in] val_span - span value.
   */
	void SetTurboRadiusOut(su2double value, unsigned short val_marker, unsigned short val_span);

  /*!
	 * \brief A total number of vertex independently from the MPI partions.
	 * \param[in] val_marker - marker value.
	 * \param[in] val_span - span value.
	 */
	unsigned long GetnTotVertexSpan(unsigned short val_marker, unsigned short val_span);

/*!
 * \brief min angular pitch independently from the MPI partions.
 * \param[in] val_marker - marker value.
 * \param[in] val_span - span value.
 */
  su2double GetMinAngularCoord(unsigned short val_marker, unsigned short val_span);

/*!
 * \brief max angular pitch independently from the MPI partions.
 * \param[in] val_marker - marker value.
 * \param[in] val_span - span value.
 */
  su2double GetMaxAngularCoord(unsigned short val_marker, unsigned short val_span);

/*!
 * \brief min Relatice angular coord independently from the MPI partions.
 * \param[in] val_marker - marker value.
 * \param[in] val_span - span value.
 */
  su2double GetMinRelAngularCoord(unsigned short val_marker, unsigned short val_span);

  /*!
	 * \brief Get the average grid velocity at a specific span for a given marker.
	 * \param[in] val_marker - marker value.
	 * \param[in] val_span - span value.
	 */
  su2double* GetAverageGridVel(unsigned short val_marker, unsigned short val_span);

};

/*! 
 * \class CMultiGridGeometry
 * \brief Class for defining the multigrid geometry, the main delicated part is the 
 *        agglomeration stage, which is done in the declaration.
 * \author F. Palacios
 */
class CMultiGridGeometry : public CGeometry {

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Level of the multigrid.
	 * \param[in] iZone - Current zone in the mesh.
	 */	
	CMultiGridGeometry(CGeometry ****geometry, CConfig **config_container, unsigned short iMesh, unsigned short iZone, unsigned short iInst);

	/*! 
	 * \brief Destructor of the class.
	 */
	~CMultiGridGeometry(void);

	/*! 
	 * \brief Determine if a CVPoint van be agglomerated, if it have the same marker point as the seed.
	 * \param[in] CVPoint - Control volume to be agglomerated.
	 * \param[in] marker_seed - Marker of the seed.
	 * \param[in] fine_grid - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \return <code>TRUE</code> or <code>FALSE</code> depending if the control volume can be agglomerated.
	 */	
	bool SetBoundAgglomeration(unsigned long CVPoint, short marker_seed, CGeometry *fine_grid, CConfig *config);

	/*! 
	 * \brief Determine if a can be agglomerated using geometrical criteria.
	 * \param[in] iPoint - Seed point.
	 * \param[in] fine_grid - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */		
	bool GeometricalCheck(unsigned long iPoint, CGeometry *fine_grid, CConfig *config);

	/*! 
	 * \brief Determine if a CVPoint van be agglomerated, if it have the same marker point as the seed.
	 * \param[in] Suitable_Indirect_Neighbors - List of Indirect Neighbours that can be agglomerated.
	 * \param[in] iPoint - Seed point.
	 * \param[in] Index_CoarseCV - Index of agglomerated point.
	 * \param[in] fine_grid - Geometrical definition of the problem.
	 */		
	void SetSuitableNeighbors(vector<unsigned long> *Suitable_Indirect_Neighbors, unsigned long iPoint, 
			unsigned long Index_CoarseCV, CGeometry *fine_grid);

	/*! 
	 * \brief Set boundary vertex.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetVertex(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Set points which surround a point.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */	
	void SetPoint_Connectivity(CGeometry *geometry);

	/*! 
	 * \brief Function declaration to avoid partially overridden classes.
	 */	
	void SetPoint_Connectivity(void);

	/*! 
	 * \brief Set the edge structure of the agglomerated control volume.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] action - Allocate or not the new elements.
	 */	
	void SetControlVolume(CConfig *config, CGeometry *geometry, unsigned short action);

	/*! 
	 * \brief Mach the near field boundary condition.
	 * \param[in] config - Definition of the particular problem.
	 */
	void MatchNearField(CConfig *config);
  
  /*!
	 * \brief Mach the near field boundary condition.
	 * \param[in] config - Definition of the particular problem.
	 */
	void MatchActuator_Disk(CConfig *config);

	/*! 
	 * \brief Mach the interface boundary condition.
	 * \param[in] config - Definition of the particular problem.
	 */
	void MatchInterface(CConfig *config);

  /*!
   * \brief Mach the periodic boundary conditions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_periodic - Index of the first periodic face in a pair.
   */
  void MatchPeriodic(CConfig *config, unsigned short val_periodic);
  
	/*! 
	 * \brief Set boundary vertex structure of the agglomerated control volume.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] action - Allocate or not the new elements.
	 */	
	void SetBoundControlVolume(CConfig *config, CGeometry *geometry, unsigned short action);

	/*! 
	 * \brief Set a representative coordinates of the agglomerated control volume.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */	
	void SetCoord(CGeometry *geometry);

        /*! 
	 * \brief Set a representative wall normal heat flux of the agglomerated control volume on a particular boundary marker.
	 * \param[in] geometry - Geometrical definition of the problem.
         * \param[in] val_marker - Index of the boundary marker.
	 */
        void SetMultiGridWallHeatFlux(CGeometry *geometry, unsigned short val_marker);

        /*! 
	 * \brief Set a representative wall temperature of the agglomerated control volume on a particular boundary marker.
	 * \param[in] geometry - Geometrical definition of the problem.
         * \param[in] val_marker - Index of the boundary marker.
	 */
        void SetMultiGridWallTemperature(CGeometry *geometry, unsigned short val_marker);

	/*!
	 * \brief Set the rotational velocity at each grid point on a coarse mesh.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] val_iZone - Index of the current zone.
	 */
	void SetRotationalVelocity(CConfig *config, unsigned short val_iZone, bool print);

/*!
 * \brief Set the rotational velocity of the points on the shroud markers to 0.0.
 * \param[in] config - Definition of the particular problem.
 */
void SetShroudVelocity(CConfig *config);

/*!
 * \brief Set the translational velocity at each grid point on a coarse mesh.
 * \param[in] config - Definition of the particular problem.
 */
void SetTranslationalVelocity(CConfig *config, unsigned short val_iZone, bool print);

	/*!
	 * \brief Set the grid velocity at each node in the coarse mesh level.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iter - Current physical time step.
	 */
	void SetGridVelocity(CConfig *config, unsigned long iter);

	/*!
	 * \brief Set the grid velocity at each node in the coarse mesh level based
	 *        on a restriction from a finer mesh.
	 * \param[in] fine_mesh - Geometry container for the finer mesh level.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetRestricted_GridVelocity(CGeometry *fine_mesh, CConfig *config);

	/*!
	 * \brief Find and store the closest neighbor to a vertex.
	 * \param[in] config - Definition of the particular problem.
	 */
	void FindNormal_Neighbor(CConfig *config);
  
	/*!
	 * \brief Indentify geometrical planes in the mesh
	 */
	void SetGeometryPlanes(CConfig *config);

	/*!
	 * \brief Get geometrical planes in the mesh
	 */
	vector<su2double> GetGeometryPlanes();

	/*!
	 * \brief Get x coords of geometrical planes in the mesh
	 */
	vector<vector<su2double> > GetXCoord();

	/*!
	 * \brief Get y coords of geometrical planes in the mesh
	 */
	vector<vector<su2double> > GetYCoord();

	/*!
	 * \brief Get z coords of geometrical planes in the mesh
	 */
	vector<vector<su2double> > GetZCoord();

	/*!
	 * \brief Get all points on a geometrical plane in the mesh
	 */
	vector<vector<unsigned long> > GetPlanarPoints();

};

/*! 
 * \struct CMultiGridQueue
 * \brief Class for a multigrid queue system
 * \author F. Palacios
 * \date Aug 12, 2012
 */
class CMultiGridQueue {
	vector<vector<unsigned long> > QueueCV; /*!< \brief Queue structure to choose the next control volume in the agglomeration process. */
	short *Priority;	/*!< \brief The priority is based on the number of pre-agglomerated neighbors. */
	bool *RightCV;	/*!< \brief In the lowest priority there are some CV that can not be agglomerated, this is the way to identify them */  
	unsigned long nPoint; /*!< \brief Total number of points. */  

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_npoint - Number of control volumes.
	 */
	CMultiGridQueue(unsigned long val_npoint);

	/*! 
	 * \brief Destructor of the class.
	 */
	~CMultiGridQueue(void);

	/*! 
	 * \brief Add a new CV to the list.
	 * \param[in] val_new_point - Index of the new point.
	 * \param[in] val_number_neighbors - Number of neighbors of the new point.
	 */
	void AddCV(unsigned long val_new_point, unsigned short val_number_neighbors);

	/*! 
	 * \brief Remove a CV from the list.
	 * \param[in] val_remove_point - Index of the control volume to be removed.
	 */
	void RemoveCV(unsigned long val_remove_point);

	/*! 
	 * \brief Change a CV from a list to a different list.
	 * \param[in] val_move_point - Index of the control volume to be moved.
	 * \param[in] val_number_neighbors - New number of neighbors of the control volume.
	 */
	void MoveCV(unsigned long val_move_point, short val_number_neighbors);

	/*! 
	 * \brief Increase the priority of the CV.
	 * \param[in] val_incr_point - Index of the control volume.
	 */
	void IncrPriorityCV(unsigned long val_incr_point);

	/*! 
	 * \brief Increase the priority of the CV.
	 * \param[in] val_red_point - Index of the control volume.
	 */
	void RedPriorityCV(unsigned long val_red_point);

	/*! 
	 * \brief Visualize the control volume queue.
	 */
	void VisualizeQueue(void);

	/*! 
	 * \brief Visualize the priority list.
	 */
	void VisualizePriority(void);

	/*! 
	 * \brief Find a new seed control volume.
	 * \return Index of the new control volume.
	 */
	long NextCV(void);

	/*! 
	 * \brief Check if the queue is empty.
	 * \return <code>TRUE</code> or <code>FALSE</code> depending if the queue is empty.
	 */
	bool EmptyQueue(void);

	/*! 
	 * \brief Total number of control volume in the queue.
	 * \return Total number of control points.
	 */
	unsigned long TotalCV(void);

	/*! 
	 * \brief Update the queue with the new control volume (remove the CV and 
	 increase the priority of the neighbors).
	 * \param[in] val_update_point - Index of the new point.
	 * \param[in] fine_grid - Fine grid geometry.
	 */
	void Update(unsigned long val_update_point, CGeometry *fine_grid);

};

#include "geometry_structure.inl"
