/*!
 * \file geometry_structure.hpp
 * \brief Headers of the main subroutines for creating the geometrical structure.
 *        The subroutines and functions are in the <i>geometry_structure.cpp</i> file.
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
  #include "cgnslib.h"
#endif
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "primal_grid_structure.hpp"
#include "dual_grid_structure.hpp"
#include "config_structure.hpp"

using namespace std;

/*! 
 * \class CGeometry
 * \brief Parent class for defining the geometry of the problem (complete geometry, 
 *        multigrid agglomerated geometry, only boundary geometry, etc..)
 * \author F. Palacios
 * \version 4.2.0 "Cardinal"
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
  unsigned long Max_GlobalPoint;  /*!< \brief Greater global point in the domain local structure. */

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
	unsigned long *nVertex;	/*!< \brief Number of vertex for each marker. */
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
#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
  idx_t * adjacency;
  idx_t * xadj;
#endif
#endif
  
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
	virtual void SetCoord_CG(void);

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
   * \param[in] mesh_filename - Name of the file where the tecplot information is going to be stored.
   * \param[in] new_file - Boolean to decide if aopen a new file or add to a old one
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetBoundSTL(char mesh_filename[MAX_STRING_SIZE], bool new_file, CConfig *config);

  
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
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetPeriodicBoundary(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] val_iZone - Index of the current zone.
	 */
	virtual void SetRotationalVelocity(CConfig *config, unsigned short val_iZone);

    /*!
     * \brief A virtual member.
     * \param[in] config - Definition of the particular problem.
     */
    virtual void SetTranslationalVelocity(CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iter - Current physical time step.
	 */
	virtual void SetGridVelocity(CConfig *config, unsigned long iter);

  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
  virtual void Set_MPI_Coord(CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
  virtual void Set_MPI_GridVel(CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
  virtual void Set_MPI_OldCoord(CConfig *config);

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
                                      su2double MinXCoord, su2double MaxXCoord, su2double *FlowVariable,
                                      vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil,
                                      vector<su2double> &Zcoord_Airfoil, vector<su2double> &Variable_Airfoil,
                                      bool original_surface, CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
  virtual su2double Compute_MaxThickness(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil, bool original_surface);
 
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
  virtual su2double Compute_AoA(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil, bool original_surface);

  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
  virtual su2double Compute_Chord(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil, bool original_surface);

  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
   * \returns The minimum value of the airfoil thickness.
	 */
	virtual su2double Compute_Thickness(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, su2double Location, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil, bool original_surface);
	
	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
   * \returns The total volume of the airfoil.
	 */
	virtual su2double Compute_Area(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil, bool original_surface);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
   * \returns The total volume of the 3D body.
	 */
  virtual su2double Compute_Volume(CConfig *config, bool original_surface);
  
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
   * \brief Update the multi-grid structure and the wall-distance.
   * \param geometry_container - Geometrical definition.
   * \param config - Config
   */
  void UpdateGeometry(CGeometry **geometry_container, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param config - Config
   */
  virtual void SetSensitivity(CConfig *config);

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
};

/*!
 * \class CPhysicalGeometry
 * \brief Class for reading a defining the primal grid which is read from the 
 *        grid file in .su2 format.
 * \author F. Palacios
 * \version 4.2.0 "Cardinal"
 */
class CPhysicalGeometry : public CGeometry {

  long *Global_to_Local_Point;				/*!< \brief Global-local indexation for the points. */
  long *Local_to_Global_Point;				/*!< \brief Local-global indexation for the points. */
  unsigned short *Local_to_Global_Marker;	/*!< \brief Local to Global marker. */
  unsigned short *Global_to_Local_Marker;	/*!< \brief Global to Local marker. */
  unsigned long *adj_counter; /*!< \brief Adjacency counter. */
  unsigned long **adjacent_elem; /*!< \brief Adjacency element list. */
  su2double* Sensitivity; /*! <\brief Vector holding the sensitivities at each point. */

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
	 * \brief Destructor of the class.
	 */
	~CPhysicalGeometry(void);
  
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
	 * \brief Set the output file for boundaries in STL CAD format
	 * \param[in] config - Definition of the particular problem.		 
	 * \param[in] mesh_filename - Name of the file where the STL 
	 *            information is going to be stored.
   * \param[in] new_file - Create a new file.
	 */
	void SetBoundSTL(char mesh_filename[MAX_STRING_SIZE], bool new_file, CConfig *config) ;

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
	 * \brief Set the rotational velocity at each node.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] val_iZone - Index of the current zone.
	 */
	void SetRotationalVelocity(CConfig *config, unsigned short val_iZone);
    
    /*!
     * \brief Set the translational velocity at each node.
     * \param[in] config - Definition of the particular problem.
     */
    void SetTranslationalVelocity(CConfig *config);

	/*! 
	 * \brief Set the grid velocity via finite differencing at each node.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetGridVelocity(CConfig *config, unsigned long iter);
  
  /*!
	 * \brief Perform the MPI communication for the grid coordinates (dynamic meshes).
	 * \param[in] config - Definition of the particular problem.
	 */
  void Set_MPI_Coord(CConfig *config);
  
  /*!
	 * \brief Perform the MPI communication for the grid velocities.
	 * \param[in] config - Definition of the particular problem.
	 */
  void Set_MPI_GridVel(CConfig *config);
  
  /*!
	 * \brief Perform the MPI communication for the grid coordinates (dynamic meshes) for restart purposes.
	 * \param[in] config - Definition of the particular problem.
	 */
  void Set_MPI_OldCoord(CConfig *config);

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
   * \brief Compute the sections of a wing.
   * \param[in] config - Definition of the particular problem.
   */
  su2double Compute_MaxThickness(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil, bool original_surface);
  
  /*!
   * \brief Compute the sections of a wing.
   * \param[in] config - Definition of the particular problem.
   */
  su2double Compute_AoA(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil, bool original_surface);
  
  /*!
   * \brief Compute the sections of a wing.
   * \param[in] config - Definition of the particular problem.
   */
  su2double Compute_Chord(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil, bool original_surface);
  
  /*!
   * \brief Find the minimum thickness of the airfoil.
   * \param[in] config - Definition of the particular problem.
   * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
   * \returns The minimum value of the airfoil thickness.
   */
  su2double Compute_Thickness(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, su2double Location, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil, bool original_surface);
  
  /*!
   * \brief Find the total volume of the airfoil.
   * \param[in] config - Definition of the particular problem.
   * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
   * \returns The total volume of the airfoil.
   */
  su2double Compute_Area(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil, bool original_surface);
  
  /*!
   * \brief Find the internal volume of the 3D body.
   * \param[in] config - Definition of the particular problem.
   * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
   * \returns The total volume of the 3D body.
   */
  su2double Compute_Volume(CConfig *config, bool original_surface);
  

  /*!
   * \brief Read the sensitivity from adjoint solution file and store it.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSensitivity(CConfig *config);

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

};

/*! 
 * \class CMultiGridGeometry
 * \brief Class for defining the multigrid geometry, the main delicated part is the 
 *        agglomeration stage, which is done in the declaration.
 * \author F. Palacios
 * \version 4.2.0 "Cardinal"
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
   * \param[in] val_iZone - Index of the current zone.
	 */
	void SetRotationalVelocity(CConfig *config, unsigned short val_iZone);
    
    /*!
     * \brief Set the translational velocity at each grid point on a coarse mesh.
     * \param[in] config - Definition of the particular problem.
     */
    void SetTranslationalVelocity(CConfig *config);

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
 * \class CPeriodicGeometry
 * \brief Class for defining a periodic boundary condition.
 * \author T. Economon, F. Palacios
 * \version 4.2.0 "Cardinal"
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
	void SetTecPlot(char config_filename[MAX_STRING_SIZE], bool new_file);

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
 * \author F. Palacios
 * \version 4.2.0 "Cardinal"
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
