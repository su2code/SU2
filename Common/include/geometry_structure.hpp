/*!
 * \file geometry_structure.hpp
 * \brief Headers of the main subroutines for creating the geometrical structure.
 *        The subroutines and functions are in the <i>geometry_structure.cpp</i> file.
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

#pragma once

#ifndef NO_MPI
#include <mpi.h>
#endif
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef NO_METIS
extern "C" {
#include "metis.h"
}
#endif

#ifndef NO_CGNS
#include "cgnslib.h"
#endif

#include "primal_grid_structure.hpp"
#include "dual_grid_structure.hpp"
#include "config_structure.hpp"

using namespace std;

/*! 
 * \class CGeometry
 * \brief Parent class for defining the geometry of the problem (complete geometry, 
 *        multigrid agglomerated geometry, only boundary geometry, etc..)
 * \author F. Palacios.
 * \version 2.0.6
 */
class CGeometry {
protected:
	unsigned long nPoint,	/*!< \brief Number of points of the mesh. */
	nPointDomain,						/*!< \brief Number of real points of the mesh. */
	nPointGhost,					/*!< \brief Number of ghost points of the mesh. */
  Global_nPoint,	/*!< \brief Total number of nodes in a simulation across all processors (including halos). */
	Global_nPointDomain,	/*!< \brief Total number of nodes in a simulation across all processors (excluding halos). */
	nElem,					/*!< \brief Number of elements of the mesh. */
  Global_nElem,	/*!< \brief Total number of elements in a simulation across all processors (all types). */
	nEdge,					/*!< \brief Number of edges of the mesh. */
	nFace,					/*!< \brief Number of faces of the mesh. */
	nElem_Storage,			/*!< \brief Storage capacity for ParaView format (domain). */
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
  nelem_wedge,          /*!< \brief Number of wedges in the mesh. */
  Global_nelem_wedge,          /*!< \brief Total number of wedges in the mesh across all processors. */
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
	bool FinestMGLevel; /*!< \brief Indicates whether the geometry class contains the finest (original) multigrid mesh. */
  unsigned long Max_GlobalPoint;  /*!< \brief Greater global point in the domain local structure. */

public:
	unsigned long *nElem_Bound_Storage;	/*!< \brief Storage capacity for ParaView format (boundaries, for each marker). */ 
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
	unsigned long *nVertex;	/*!< \brief Number of vertex for each marker. */
	unsigned short nCommLevel;		/*!< \brief Number of non-blocking communication levels. */
	vector<unsigned long> PeriodicPoint[MAX_NUMBER_PERIODIC][2];			/*!< \brief PeriodicPoint[Periodic bc] and return the point that 
																			 must be sent [0], and the image point in the periodic bc[1]. */
	vector<unsigned long> PeriodicElem[MAX_NUMBER_PERIODIC];				/*!< \brief PeriodicElem[Periodic bc] and return the elements that 
																			 must be sent. */
	vector<unsigned long> OldBoundaryElems[MAX_NUMBER_MARKER];  /*!< \brief Vector of old boundary elements. */

  
  
  vector<unsigned long> SendTransfLocal[MAX_NUMBER_DOMAIN];	/*!< \brief Vector to store the type of transformation for this
                                                             send point. */
  vector<unsigned long> ReceivedTransfLocal[MAX_NUMBER_DOMAIN];	/*!< \brief Vector to store the type of transformation for this
                                                                 received point. */
	vector<unsigned long> SendDomainLocal[MAX_NUMBER_DOMAIN];								/*!< \brief SendDomain[from domain][to domain] and return the
                                                                           point index of the node that must me sended. */
	vector<unsigned long> ReceivedDomainLocal[MAX_NUMBER_DOMAIN];								/*!< \brief SendDomain[from domain][to domain] and return the
                                                                               point index of the node that must me sended. */
  short *Marker_All_SendRecv;
  
	/*--- Create vectors and distribute the values among the different planes queues ---*/
	vector<vector<double> > Xcoord_plane; /*!< \brief Vector containing x coordinates of new points appearing on a single plane */
	vector<vector<double> > Ycoord_plane; /*!< \brief Vector containing y coordinates of  new points appearing on a single plane */
	vector<vector<double> > Zcoord_plane; 	/*!< \brief Vector containing z coordinates of  new points appearing on a single plane */
	vector<vector<double> > FaceArea_plane; /*!< \brief Vector containing area/volume associated with  new points appearing on a single plane */
	vector<vector<unsigned long> > Plane_points; /*!< \brief Vector containing points appearing on a single plane */

	vector<double> XCoordList;	/*!< \brief Vector containing points appearing on a single plane */
	CPrimalGrid*** newBound;            /*!< \brief Boundary vector for new periodic elements (primal grid information). */
	unsigned long *nNewElem_Bound;			/*!< \brief Number of new periodic elements of the boundary. */

	/*! 
	 * \brief Constructor of the class.
	 */
	CGeometry(void);

	/*! 
	 * \brief Destructor of the class.
	 */
	virtual ~CGeometry(void);

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
	unsigned long GetnVertex(unsigned short val_marker);

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
	double Point2Plane_Distance(double *Coord, double *iCoord, double *jCoord, double *kCoord);

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
	 * \param[in] val_ndim - Number of dimensions.
	 */
	void SetnDim(unsigned short val_ndim);

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
	 * \brief Set the number of storage for boundary elements.
	 * \param[in] val_marker - Marker of the boundary.
	 * \param[in] val_nelem_bound - Number of boundary elements.
	 */	
	void SetnElem_Bound_Storage(unsigned short val_marker, unsigned long val_nelem_bound);

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
	 * \brief Get the number of storage boundary elements.
	 * \param[in] val_marker - Marker of the boundary.
	 */
	unsigned long GetnElem_Bound_Storage(unsigned short val_marker);

	/*! 
	 * \brief Set the number of elements in vtk fortmat.
	 * \param[in] val_nelem_storage - Number of elements
	 */
	void SetnElem_Storage(unsigned long val_nelem_storage);

	/*! 
	 * \brief Get the number of elements in vtk fortmat.
	 */	
	unsigned long GetnElem_Storage(void);

  /*!
	 * \brief Get the number of elements in vtk fortmat.
	 */
	unsigned long GetMax_GlobalPoint(void);
  
	/*!
	 * \brief Get boolean for whether this is the finest (original) multigrid mesh level.
	 * \return <code>TRUE</code> if this is the finest multigrid mesh level; otherwise <code>FALSE</code>.
	 */
	bool GetFinestMGLevel(void);

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
	virtual void SetWall_Distance(CConfig *config);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.		 
	 */
	virtual void SetPositive_ZArea(CConfig *config);

	/*! 
	 * \brief A virtual member.
	 */
	virtual void SetEsuP(void);

	/*! 
	 * \brief A virtual member.
	 */	
	virtual void SetPsuP(void);

	/*! 
	 * \brief A virtual member.
	 */		
	virtual void SetEsuE(void);

	/*! 
	 * \brief A virtual member.
	 */
	void SetEdges(void);

	/*! 
	 * \brief A virtual member.
	 */
	void SetFaces(void);

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
	 */
	virtual void SetVertex(void);

	/*! 
	 * \brief A virtual member.
	 */		
	virtual void SetCG(void);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] action - Allocate or not the new elements.		 
	 */
	virtual void SetControlVolume(CConfig *config, unsigned short action);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void MatchNearField(CConfig *config);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void MatchInterface(CConfig *config);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] geometry_donor - Geometry of the donor zone.
	 * \param[in] config_donor - Definition of the donor problem.
	 */
	virtual void MatchZone(CConfig *config, CGeometry *geometry_donor, CConfig *config_donor, 
			unsigned short val_iZone, unsigned short val_nZone);

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
	virtual void SetTecPlot(char config_filename[200]);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.		 
	 * \param[in] mesh_filename - Name of the file where the tecplot information is going to be stored.
	 */
	virtual void SetBoundTecPlot(CConfig *config, char mesh_filename[200]);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.		 
	 */
	virtual void Check_Orientation(CConfig *config);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.		 
	 */
	virtual void SetColorGrid(CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
  virtual void DivideConnectivity(CConfig *config, unsigned short Elem_Type);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.		 
	 */
	virtual void SetPeriodicBoundary(CConfig *config);

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
	 */	
	virtual void SetCoord(CGeometry *geometry);

	/*! 
	 * \brief A virtual member.
	 * \param[in] val_nSmooth - Number of smoothing iterations.
	 * \param[in] val_smooth_coeff - Relaxation factor.
	 * \param[in] config - Definition of the particular problem.
	 */	
	virtual void SetCoord_Smoothing(unsigned short val_nSmooth, double val_smooth_coeff, CConfig *config);

	/*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */	
	virtual void SetPsuP(CGeometry *geometry);

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
	 * \param[in] val_mesh_out_filename - Name of the output file.
	 */
	virtual void SetMeshFile(CConfig *config, string val_mesh_out_filename, string val_mesh_in_filename);
  
	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] mesh_vtk - Name of the vtk file.
	 * \param[in] mesh_su2 - Name of the su2 file.
	 * \param[in] nslices - Number of slices of the 2D configuration.
	 */	
	virtual void Set3D_to_2D(CConfig *config, char mesh_vtk[200], char mesh_su2[200], unsigned short nslices);

	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetBoundSensitivity(CConfig *config);

	/*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetPeriodicBoundary(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetRotationalVelocity(CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iter - Current physical time step.
	 */
	virtual void SetGridVelocity(CConfig *config, unsigned long iter);

	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iter - Current physical time step.
	 */
	virtual void SetRestricted_GridVelocity(CGeometry *fine_mesh, CConfig *config, unsigned long iter);

	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void FindSharpEdges(CConfig *config);

  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
   * \returns The maximum value of the airfoil thickness.
	 */
	virtual double GetMaxThickness(CConfig *config, bool original_surface);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
   * \returns The minimum value of the airfoil thickness.
	 */
	virtual double GetMinThickness(CConfig *config, bool original_surface);
	
	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
   * \returns The total volume of the airfoil.
	 */
	virtual double GetTotalVolume(CConfig *config, bool original_surface);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
   * \returns The clearance height of the airfoil.
	 */
	virtual double GetClearance(CConfig *config, bool original_surface);
	
	/*! 
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void FindNormal_Neighbor(CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_ipoint - Global point.
	 * \returns Local index that correspond with the global index.
	 */
	virtual long GetGlobal_to_Local_Point(long val_ipoint);

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
	 * \returns Total number of elements in a simulation across all processors.
	 */
	virtual unsigned long GetGlobal_nElem();
  
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
	 * \returns Total number of wedge elements in a simulation across all processors.
	 */
	virtual unsigned long GetGlobal_nElemWedg();
  
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
	 * \return Number of wedge elements.
	 */
	virtual unsigned long GetnElemWedg();
  
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
	virtual vector<double> GetGeometryPlanes();

	/*!
	 * \brief Get x coords of geometrical planes in the mesh
	 */
	virtual vector<vector<double> > GetXCoord();

	/*!
	 * \brief Get y coords of geometrical planes in the mesh
	 */
	virtual vector<vector<double> > GetYCoord();

	/*!
	 * \brief Get z coords of geometrical planes in the mesh
	 */
	virtual vector<vector<double> > GetZCoord();

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
	void SetSpline(vector<double> &x, vector<double> &y, unsigned long n, double yp1, double ypn, vector<double> &y2);
	
	/*!
	 * \brief Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai’s in order), 
	          and given the array y2a[1..n], which is the output from spline above, and given a value of 
	          x, this routine returns a cubic-spline interpolated value y.
         	  Numerical Recipes: The Art of Scientific Computing, Third Edition in C++.
	 * \returns The interpolated value of for x.
	 */
	double GetSpline(vector<double> &xa, vector<double> &ya, vector<double> &y2a, unsigned long n, double x);
	
};

/*! 
 * \class CPhysicalGeometry
 * \brief Class for reading a defining the primal grid which is read from the 
 *        grid file in .su2 format.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CPhysicalGeometry : public CGeometry {

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
	CPhysicalGeometry(CConfig *config, string val_mesh_filename, unsigned short val_format, 
			unsigned short val_iZone, unsigned short val_nZone);

	/*! 
	 * \brief Destructor of the class.
	 */
	~CPhysicalGeometry(void);
  
  /*!
	 * \brief Reads the geometry of the grid and adjust the boundary
	 *        conditions with the configuration file.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_mesh_filename - Name of the file with the grid information.
	 * \param[in] val_format - Format of the file with the grid information.
	 * \param[in] val_iZone - Domain to be read from the grid file.
	 * \param[in] val_nZone - Total number of domains in the grid file.
	 */
	void SU2_Format(CConfig *config, string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone);
  
  /*!
	 * \brief Reads the geometry of the grid and adjust the boundary
	 *        conditions with the configuration file.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_mesh_filename - Name of the file with the grid information.
	 * \param[in] val_format - Format of the file with the grid information.
	 * \param[in] val_iZone - Domain to be read from the grid file.
	 * \param[in] val_nZone - Total number of domains in the grid file.
	 */
	void CGNS_Format(CConfig *config, string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone);
  
  /*!
	 * \brief Reads the geometry of the grid and adjust the boundary
	 *        conditions with the configuration file.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_mesh_filename - Name of the file with the grid information.
	 * \param[in] val_format - Format of the file with the grid information.
	 * \param[in] val_iZone - Domain to be read from the grid file.
	 * \param[in] val_nZone - Total number of domains in the grid file.
	 */
	void NETCDF_Format(CConfig *config, string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone);

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
	 * \brief Set elements which surround a point.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetWall_Distance(CConfig *config);

	/*! 
	 * \brief Compute positive z area for adimensionalization.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPositive_ZArea(CConfig *config);

	/*! 
	 * \brief Set elements which surround a point.
	 */
	void SetEsuP(void);

	/*! 
	 * \brief Set points which surround a point.
	 */
	void SetPsuP(void);
  
  /*!
	 * \brief Set points which surround a point. Special version for FEA
   * grid deformation that performs element divisions internally.
	 */
	void SetPsuP_FEA(void);

	/*!
	 * \brief Function declaration to avoid partially overridden classes.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */
	void SetPsuP(CGeometry *geometry);

	/*! 
	 * \brief Set elements which surround an element.
	 */
	void SetEsuE(void);

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
	 * \brief Set the center of gravity of the face, elements and edges.
	 */
	void SetCG(void);

	/*! 
	 * \brief Set the edge structure of the control volume.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] action - Allocate or not the new elements.
	 */
	void SetControlVolume(CConfig *config, unsigned short action);

	/*! 
	 * \brief Mach the near field boundary condition.
	 * \param[in] config - Definition of the particular problem.
	 */
	void MatchNearField(CConfig *config);

	/*! 
	 * \brief Mach the interface boundary condition.
	 * \param[in] config - Definition of the particular problem.
	 */
	void MatchInterface(CConfig *config);

	/*! 
	 * \brief Mach the interface boundary condition.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] geometry_donor - Geometry of the donor zone.
	 * \param[in] config_donor - Definition of the donor problem.
	 */
	void MatchZone(CConfig *config, CGeometry *geometry_donor, CConfig *config_donor, 
			unsigned short val_iZone, unsigned short val_nZone);

	/*! 
	 * \brief Set boundary vertex structure of the control volume.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] action - Allocate or not the new elements.
	 */
	void SetBoundControlVolume(CConfig *config, unsigned short action);

	/*! 
	 * \brief Set the Tecplot file.
	 * \param[in] config_filename - Name of the file where the Tecplot 
	 *            information is going to be stored.
	 */
	void SetTecPlot(char config_filename[200]);

	/*! 
	 * \brief Set the output file for boundaries in Tecplot
	 * \param[in] config - Definition of the particular problem.		 
	 * \param[in] mesh_filename - Name of the file where the Tecplot 
	 *            information is going to be stored.
	 */
	void SetBoundTecPlot(CConfig *config, char mesh_filename[200]);

	/*! 
	 * \brief Set the output file for boundaries in STL CAD format
	 * \param[in] config - Definition of the particular problem.		 
	 * \param[in] mesh_filename - Name of the file where the STL 
	 *            information is going to be stored.
	 */
	void SetBoundSTL(CConfig *config, char mesh_filename[200]);

	/*! 
	 * \brief Check the volume element orientation.
	 * \param[in] config - Definition of the particular problem.		 
	 */
	void Check_Orientation(CConfig *config);

	/*! 
	 * \brief Set the domains for grid grid partitioning.
	 * \param[in] config - Definition of the particular problem.		 
	 */
	void SetColorGrid(CConfig *config);
  
	/*!
	 * \brief Set the rotational velocity at each grid point.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetRotationalVelocity(CConfig *config);

	/*! MC - fill this in - 7/11/12
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetGridVelocity(CConfig *config, unsigned long iter);

	/*! 
	 * \brief Set the periodic boundary conditions.
	 * \param[in] config - Definition of the particular problem.		 
	 */
	void SetPeriodicBoundary(CConfig *config);

	/*! 
	 * \brief Do an implicit smoothing of the grid coordinates.
	 * \param[in] val_nSmooth - Number of smoothing iterations.
	 * \param[in] val_smooth_coeff - Relaxation factor.
	 * \param[in] config - Definition of the particular problem.		 
	 */	
	void SetCoord_Smoothing(unsigned short val_nSmooth, double val_smooth_coeff, CConfig *config);

	/*! 
	 * \brief Write the .su2 file.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_mesh_out_filename - Name of the output file.
	 */	
	void SetMeshFile(CConfig *config, string val_mesh_out_filename);
  
  /*!
	 * \brief Write the .su2 file, with new domain coordinates
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_mesh_out_filename - Name of the output file.
	 */
	void SetMeshFile(CConfig *config, string val_mesh_out_filename, string val_mesh_in_filename);

	/*! 
	 * \brief Create a 2D mesh using a 3D mesh with symmetries.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] mesh_vtk - Name of the vtk file.
	 * \param[in] mesh_su2 - Name of the su2 file.
	 * \param[in] nslices - Number of slices of the 2D configuration.
	 */	
	void Set3D_to_2D(CConfig *config, char mesh_vtk[200], char mesh_su2[200], unsigned short nslices);

	/*! 
	 * \brief Compute some parameters about the grid quality.
	 * \param[out] statistics - Information about the grid quality, statistics[0] = (r/R)_min, statistics[1] = (r/R)_ave.		 
	 */	
	void GetQualityStatistics(double *statistics);

	/*!
	 * \brief Find and store all vertices on a sharp corner in the geometry.
	 * \param[in] config - Definition of the particular problem.
	 */
	void FindSharpEdges(CConfig *config);

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
	 * \brief Retrieve total number of wedge elements in a simulation across all processors.
	 * \returns Total number of wedge elements in a simulation across all processors.
	 */
	unsigned long GetGlobal_nElemWedg();
  
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
	 * \brief Get number of wedge elements.
	 * \return Number of wedge elements.
	 */
	unsigned long GetnElemWedg();
  
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
	vector<double> GetGeometryPlanes();

	/*!
	 * \brief Get x coords of geometrical planes in the mesh
	 */
	vector<vector<double> > GetXCoord();

	/*!
	 * \brief Get y coords of geometrical planes in the mesh
	 */
	vector<vector<double> > GetYCoord();

	/*!
	 * \brief Get z coords of geometrical planes in the mesh
	 */
	vector<vector<double> > GetZCoord();

	/*!
	 * \brief Get all points on a geometrical plane in the mesh
	 */
	vector<vector<unsigned long> > GetPlanarPoints();
};

/*! 
 * \class CMultiGridGeometry
 * \brief Class for defining the multigrid geometry, the main delicated part is the 
 *        agglomeration stage, which is done in the declaration.
 * \author F. Palacios.
 * \version 2.0.6
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
	CMultiGridGeometry(CGeometry ***geometry, CConfig **config_container, unsigned short iMesh, unsigned short iZone);

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
	void SetPsuP(CGeometry *geometry);

	/*! 
	 * \brief Function declaration to avoid partially overridden classes.
	 */	
	void SetPsuP(void);

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
	 * \brief Mach the interface boundary condition.
	 * \param[in] config - Definition of the particular problem.
	 */
	void MatchInterface(CConfig *config);

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
	 * \brief Set the rotational velocity at each grid point on a coarse mesh.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetRotationalVelocity(CConfig *config);

	/*!
	 * \brief Set the grid velocity at each node in the coarse mesh level.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iter - Current physical time step.
	 */
	void SetGridVelocity(CConfig *config, unsigned long iter);

	/*!
	 * \brief Set the grid velocity at each node in the coarse mesh level based
	 *        on a restriction from a finer mesh (needed for the unsteady adjoint).
	 * \param[in] fine_mesh - Geometry container for the finer mesh level.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iter - Current physical time step.
	 */
	void SetRestricted_GridVelocity(CGeometry *fine_mesh, CConfig *config, unsigned long iter);

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
	vector<double> GetGeometryPlanes();

	/*!
	 * \brief Get x coords of geometrical planes in the mesh
	 */
	vector<vector<double> > GetXCoord();

	/*!
	 * \brief Get y coords of geometrical planes in the mesh
	 */
	vector<vector<double> > GetYCoord();

	/*!
	 * \brief Get z coords of geometrical planes in the mesh
	 */
	vector<vector<double> > GetZCoord();

	/*!
	 * \brief Get all points on a geometrical plane in the mesh
	 */
	vector<vector<unsigned long> > GetPlanarPoints();

};

/*! 
 * \class CBoundaryGeometry
 * \brief Class for only defining the boundary of the geometry, this class is only 
 *        used in case we are not interested in the volumetric grid.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CBoundaryGeometry : public CGeometry {
  
public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_mesh_filename - Name of the file with the grid information, be careful 
	 *            because as input file we don't use a .su2, in this case we use a .csv file.
	 * \param[in] val_format - Format of the file with the grid information.
	 */
	CBoundaryGeometry(CConfig *config, string val_mesh_filename, unsigned short val_format);

	/*! 
	 * \brief Destructor of the class.
	 */
	~CBoundaryGeometry(void);

	/*! 
	 * \brief Set boundary vertex.
	 */
	void SetVertex(void);

	/*! 
	 * \brief Compute the boundary geometrical structure.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] action - Allocate or not the new elements.
	 */
	void SetBoundControlVolume(CConfig *config, unsigned short action);

	/*! 
	 * \brief Read the sensitivity from an input file.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetBoundSensitivity(CConfig *config);
	
	/*! 
	 * \brief Find the maximum thickness of the airfoil.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
   * \returns The maximum value of the airfoil thickness.
	 */
  double GetMaxThickness(CConfig *config, bool original_surface);
	
  /*!
	 * \brief Find the minimum thickness of the airfoil.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
   * \returns The minimum value of the airfoil thickness.
	 */
  double GetMinThickness(CConfig *config, bool original_surface);
  
	/*! 
	 * \brief Find the total volume of the airfoil.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
   * \returns The total volume of the airfoil.
	 */
  double GetTotalVolume(CConfig *config, bool original_surface);
  
  /*!
	 * \brief Find the clearance height of the airfoil.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
   * \returns The clearance height of the airfoil.
	 */
  double GetClearance(CConfig *config, bool original_surface);
	
};

/*! 
 * \class CDomainGeometry
 * \brief Class for defining an especial kind of grid used in the partioning stage.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CDomainGeometry : public CGeometry {
	long *Global_to_Local_Point;				/*!< \brief Global-local indexation for the points. */
	unsigned long *Local_to_Global_Point;				/*!< \brief Local-global indexation for the points. */
	unsigned short *Local_to_Global_Marker;	/*!< \brief Local to Global marker. */
	unsigned short *Global_to_Local_Marker;	/*!< \brief Global to Local marker. */

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_domain - Number of domains for parallelization purposes.	 
	 */
	CDomainGeometry(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class.
	 */
	~CDomainGeometry(void);
  
  /*!
	 * \brief Constructor of the class.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_domain - Number of domains for parallelization purposes.
	 */
	void SetDomainSerial(CGeometry *geometry, CConfig *config, unsigned short val_domain);
  
	/*! 
	 * \brief Set the send receive boundaries of the grid.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_domain - Number of domains for parallelization purposes.	 
	 */
	void SetSendReceive(CConfig *config);

	/*! 
	 * \brief Set the Tecplot file.
	 * \param[in] config_filename - Name of the file where the Tecplot
	 *            information is going to be stored.
	 */
	void SetTecPlot(char config_filename[200]);

	/*! 
	 * \brief Write the .su2 file.
	 * \param[in] config - Definition of the particular problem.		 
	 * \param[in] val_mesh_out_filename - Name of the output file.
	 */
	void SetMeshFile(CConfig *config, string val_mesh_out_filename);

	/*!
	 * \brief Get the local index that correspond with the global numbering index.
	 * \param[in] val_ipoint - Global point.
	 * \returns Local index that correspond with the global index.
	 */
	long GetGlobal_to_Local_Point(long val_ipoint);

	/*!
	 * \brief Get the local marker that correspond with the global marker.
	 * \param[in] val_ipoint - Global marker.
	 * \returns Local marker that correspond with the global index.
	 */
	unsigned short GetGlobal_to_Local_Marker(unsigned short val_imarker);

};

/*! 
 * \class CPeriodicGeometry
 * \brief Class for defining a periodic boundary condition.
 * \author T. Economon, F. Palacios.
 * \version 2.0.6
 */
class CPeriodicGeometry : public CGeometry {
	CPrimalGrid*** newBoundPer;            /*!< \brief Boundary vector for new periodic elements (primal grid information). */
	unsigned long *nNewElem_BoundPer;			/*!< \brief Number of new periodic elements of the boundary. */

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CPeriodicGeometry(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class.
	 */
	~CPeriodicGeometry(void);

	/*! 
	 * \brief Set the periodic boundaries of the grid.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPeriodicBoundary(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Set the Tecplot file.
	 * \param[in] config_filename - Name of the file where the Tecplot 
	 *            information is going to be stored.
	 */
	void SetTecPlot(char config_filename[200]);

	/*! 
	 * \brief Write the .su2 file.
	 * \param[in] config - Definition of the particular problem.		 
	 * \param[in] val_mesh_out_filename - Name of the output file.
	 */
	void SetMeshFile(CGeometry *geometry, CConfig *config, string val_mesh_out_filename);
};

/*! 
 * \struct CMultiGridQueue
 * \brief Class for a multigrid queue system
 * \author F. Palacios.
 * \version 2.0.6
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
