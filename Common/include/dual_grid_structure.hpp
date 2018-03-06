/*!
 * \file dual_grid_structure.hpp
 * \brief Headers of the main subroutines for doing the complete dual grid structure.
 *        The subroutines and functions are in the <i>dual_grid_structure.cpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 6.0.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>

#include "config_structure.hpp"

using namespace std;

/*! 
 * \class CDualGrid
 * \brief Class for controlling the dual volume definition. The dual volume is compose by 
 *        three main elements: points, edges, and vertices.
 * \author F. Palacios
 */
class CDualGrid{
protected:
	static unsigned short nDim; /*!< \brief Number of dimensions of the problem. */
	
public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.		 
	 */
	CDualGrid(unsigned short val_nDim);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	virtual ~CDualGrid(void);
	
	/*! 
	 * \brief A pure virtual member.
	 */
	virtual su2double *GetCoord(void) = 0;
	
	/*! 
	 * \brief A pure virtual member.
	 * \param[in] val_coord - Coordinate of the point.		 
	 */
	virtual void SetCoord(su2double *val_coord) = 0;
	
	/*! 
	 * \brief A pure virtual member.
	 * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
	 * \param[in] val_coord_FaceElem_CG - Coordinates of the centre of gravity of the face of an element.
	 * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG) = 0;
	
	/*! 
	 * \overload
	 * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
	 * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG) = 0;
	
	/*! 
	 * \brief A pure virtual member.
	 * \param[in] val_normal - Coordinates of the normal.
	 */
	virtual void GetNormal(su2double *val_normal) = 0;
	
	/*! 
	 * \brief A pure virtual member.
	 */
	virtual su2double *GetNormal(void) = 0;
	
	/*! 
	 * \brief A pure virtual member.
	 * \param[in] val_face_normal - Coordinates of the normal.
	 */	
	virtual void SetNormal(su2double *val_face_normal) = 0;

	/*! 
	 * \brief A pure virtual member.
	 */
	virtual unsigned short GetnNodes(void) = 0;

	/*! 
	 * \brief A pure virtual member.
	 */	
	virtual void SetZeroValues(void) = 0;
	
	/*! 
	 * \brief A pure virtual member.
	 * \param[in] val_face_normal - Normal vector to be added.
	 */	
	virtual void AddNormal(su2double *val_face_normal) = 0;
};

/*! 
 * \class CPoint
 * \brief Class for point definition (including control volume definition).
 * \author F. Palacios
 */
class CPoint : public CDualGrid {
private:
  unsigned short nElem,               /*!< \brief Number of elements that set up the control volume. */
  nPoint;                             /*!< \brief Number of points that set up the control volume  */
  vector<long> Elem;                  /*!< \brief Elements that set up a control volume around a node. */
  vector<unsigned long> Point;        /*!< \brief Points surrounding the central node of the control volume. */
  vector<long> Edge;                  /*!< \brief Edges that set up a control volume. */
  su2double *Volume;                  /*!< \brief Volume or Area of the control volume in 3D and 2D. */
  bool Domain,                        /*!< \brief Indicates if a point must be computed or belong to another boundary */
  Boundary,                           /*!< \brief To see if a point belong to the boundary (including MPI). */
  PhysicalBoundary,                   /*!< \brief To see if a point belong to the physical boundary (without includin MPI). */
  SolidBoundary;                      /*!< \brief To see if a point belong to the physical boundary (without includin MPI). */
  long *Vertex;                       /*!< \brief Index of the vertex that correspond which the control volume (we need one for each marker in the same node). */
  su2double *Coord,                   /*!< \brief vector with the coordinates of the node. */
  *Coord_Old,                         /*!< \brief Old coordinates vector for geometry smoothing. */
  *Coord_Sum,                         /*!< \brief Sum of coordinates vector for geometry smoothing. */
  *Coord_n,                           /*!< \brief Coordinates at time n for use with dynamic meshes. */
  *Coord_n1,                          /*!< \brief Coordinates at time n-1 for use with dynamic meshes. */
  *Coord_p1;                          /*!< \brief Coordinates at time n+1 for use with dynamic meshes. */
  su2double *GridVel;                 /*!< \brief Velocity of the grid for dynamic mesh cases. */
  su2double **GridVel_Grad;           /*!< \brief Gradient of the grid velocity for dynamic meshes. */
  unsigned long Parent_CV;            /*!< \brief Index of the parent control volume in the agglomeration process. */
  unsigned short nChildren_CV;        /*!< \brief Number of children in the agglomeration process. */
  vector<unsigned long> Children_CV;  /*!< \brief Index of the children control volumes in the agglomeration process. */
  bool Agglomerate_Indirect,					/*!< \brief This flag indicates if the indirect points can be agglomerated. */
  Agglomerate;                        /*!< \brief This flag indicates if the element has been agglomerated. */
  bool Move;                          /*!< \brief This flag indicates if the point is going to be move in the grid deformation process. */
  unsigned short color;               /*!< \brief Color of the point in the partitioning strategy. */
  su2double Wall_Distance;            /*!< \brief Distance to the nearest wall. */
  su2double SharpEdge_Distance;       /*!< \brief Distance to a sharp edge. */
  su2double Curvature;                /*!< \brief Value of the surface curvature (SU2_GEO). */
  unsigned long GlobalIndex;          /*!< \brief Global index in the parallel simulation. */
  unsigned short nNeighbor;           /*!< \brief Number of neighbors. */
  bool Flip_Orientation;              /*!< \brief Flip the orientation of the normal. */

public:
	
	/*! 
	 * \brief Constructor of the class. 
	 * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_globalindex Global index in the parallel simulation.
	 * \param[in] config - Definition of the particular problem.
	 */
	CPoint(unsigned short val_nDim, unsigned long val_globalindex, CConfig *config);
	
	/*! 
	 * \overload
	 * \param[in] val_coord_0 First coordinate of the point.		 
	 * \param[in] val_coord_1 Second coordinate of the point.
	 * \param[in] val_globalindex Global index in the parallel simulation.
	 * \param[in] config - Definition of the particular problem.
	 */
	CPoint(su2double val_coord_0, su2double val_coord_1, unsigned long val_globalindex, CConfig *config);
	
	/*! 
	 * \overload
	 * \param[in] val_coord_0 First coordinate of the point.		 
	 * \param[in] val_coord_1 Second coordinate of the point.		
	 * \param[in] val_coord_2 Third coordinate of the point.
	 * \param[in] val_globalindex Global index in the parallel simulation.
	 * \param[in] config - Definition of the particular problem.
	 */
	CPoint(su2double val_coord_0, su2double val_coord_1, su2double val_coord_2, unsigned long val_globalindex, CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CPoint(void);
	
	/*! 
	 * \brief For parallel computation, its indicates if a point must be computed or not.
	 * \param[in] val_domain - <code>TRUE</code> if the point belong to the domain; otherwise <code>FALSE</code>.
	 */
	void SetDomain(bool val_domain);	
	
	/*! 
	 * \brief For parallel computation, its indicates if a point must be computed or not.
	 * \return <code>TRUE</code> if the node belong to the physical domain; otherwise <code>FALSE</code>.
	 */
	bool GetDomain(void);
	
	/*! 
	 * \brief Set the value of the distance to the nearest wall.
	 * \param[in] val_distance - Value of the distance.
	 */
	void SetWall_Distance(su2double val_distance);
  
  /*!
	 * \brief Set the value of the distance to a sharp edge.
	 * \param[in] val_distance - Value of the distance.
	 */
	void SetSharpEdge_Distance(su2double val_distance);
	
	/*! 
	 * \brief Get the value of the distance to the nearest wall.
	 * \return Value of the distance to the nearest wall.
	 */
	su2double GetWall_Distance(void);
	
  /*!
	 * \brief Set the value of the curvature at a surface node.
	 * \param[in] val_curvature - Value of the curvature.
	 */
	void SetCurvature(su2double val_curvature);
	
	/*!
	 * \brief Get the value of the curvature at a surface node.
	 * \return Value of the curvature.
	 */
	su2double GetCurvature(void);
  
  /*!
	 * \brief Get the value of the distance to a sharp edge
	 * \return Value of the distance to the nearest wall.
	 */
	su2double GetSharpEdge_Distance(void);
  
	/*! 
	 * \brief Set the number of elements that compose the control volume.
	 * \param[in] val_nElem - Number of elements that make the control volume around a node.
	 */
	void SetnElem(unsigned short val_nElem);
	
	/*! 
	 * \brief Set the number of points that compose the control volume.
	 * \param[in] val_nPoint - Number of points that compose the control volume (points surrounding points).
	 */
	void SetnPoint(unsigned short val_nPoint);
	
	/*! 
	 * \brief Get the coordinates dor the control volume.
	 * \param[in] val_dim - Number of dimensions of the problem.
	 * \return Coordinate that correspond with <i>val_dim</i>.
	 */
	su2double GetCoord(unsigned short val_dim);
	
	/*! 
	 * \brief Get the coordinates of the control volume.
	 * \return pointer to the coordinate of the point.
	 */
	su2double *GetCoord(void);
	
	/*! 
	 * \brief Set the coordinates for the control volume.
	 * \param[in] val_dim - Position to store the coordinate.		 
	 * \param[in] val_coord - Coordinate for val_dim.			 
	 */
	void SetCoord(unsigned short val_dim, su2double val_coord);
  
  /*!
	 * \brief Get the coordinates of the control volume.
	 * \return pointer to the coordinate of the point.
	 */
	bool GetFlip_Orientation(void);
	
	/*!
	 * \brief Set the coordinates for the control volume.
	 * \param[in] val_dim - Position to store the coordinate.
	 * \param[in] val_coord - Coordinate for val_dim.
	 */
	void SetFlip_Orientation(void);
  
  /*!
	 * \brief Set the coordinates for the control volume.
	 * \param[in] val_dim - Position to store the coordinate.
	 * \param[in] val_coord - Coordinate for val_dim.
	 */
	void AddCoord(unsigned short val_dim, su2double val_coord);
  
	/*! 
	 * \overload
	 * \param[in] val_coord - Coordinate of the point.		 
	 */
	void SetCoord(su2double *val_coord);
	
	/*! 
	 * \brief Get the number of elements that compose the control volume.
	 * \return Number of elements that compose the control volume.
	 */
	unsigned short GetnElem(void);
	
	/*! 
	 * \brief Get the number of points that compose the control volume.
	 * \return Number of points that compose the control volume.
	 */
	unsigned short GetnPoint(void);
	
	/*! 
	 * \brief Set the elements that set the control volume.
	 * \param[in] val_elem - Element to be added.		 
	 */
	void SetElem(unsigned long val_elem);
  
  /*!
	 * \brief Reset the elements of a control volume.
	 */
	void ResetElem(void);
  
  /*!
	 * \brief Reset the points that compose the control volume.
	 */
	void ResetPoint(void);

	/*! 
	 * \brief Set the points that compose the control volume.
	 * \param[in] val_point - Point to be added.		 
	 */
	void SetPoint(unsigned long val_point);
	
	/*! 
	 * \brief Set the edges that compose the control volume.
	 * \param[in] val_edge - Edge to be added.
	 * \param[in] val_nEdge - Position in which is going to be stored the edge for each control volume.
	 */
	void SetEdge(long val_edge, unsigned short val_nEdge);
	
	/*! 
	 * \brief Set the boundary vertex that compose the control volume.
	 * \param[in] val_vertex - Vertex to be added.
	 * \param[in] val_nMarker - Marker of the vertex to be added (position where is going to be stored).
	 */
	void SetVertex(long val_vertex, unsigned short val_nMarker);
	
	/*! 
	 * \brief Get all the elements that compose the control volume.
	 * \param[in] val_elem - Position where the element is stored.	 
	 * \return Index of the element.
	 */
	unsigned long GetElem(unsigned short val_elem);
	
	/*! 
	 * \brief Get all the points that compose the control volume.
	 * \param[in] val_point - Position where the point is stored.
	 * \return Index of the point.
	 */
	unsigned long GetPoint(unsigned short val_point);
	
	/*! 
	 * \brief Get all the edges that compose the control volume.
	 * \param[in] val_edge - Position where the edge is stored.
	 * \return Index of the edge.
	 */
	long GetEdge(unsigned short val_edge);
	
	/*! 
	 * \brief Get the vertex that compose the control volume for a marker.
	 * \param[in] val_marker - Position where the vertex is stored.		 
	 * \return Index of the vertex.
	 */
	long GetVertex(unsigned short val_marker);
	
	/*! 
	 * \brief Adds some area or volume of the CV.
	 * \param[in] val_Volume - Local volume to be added to the total one.
	 */
	void AddVolume(su2double val_Volume);
	
	/*! 
	 * \brief Get area or volume of the control volume.
	 * \return Area or volume of the control volume.
	 */
	su2double GetVolume(void);
	
	/*! 
	 * \brief Get information about the movement of the node.
	 * \return <code>TRUE</code> if the point is going to be moved; otherwise <code>FALSE</code>.
	 */
	bool GetMove(void);
	
	/*! 
	 * \brief Set if a point belong to the boundary.
	 * \note It also create the structure to store the vertex.
	 * \param[in] val_nmarker - Max number of marker.
	 */		
	void SetBoundary(unsigned short val_nmarker);
  
  /*!
	 * \brief Reset the boundary of a control volume.
	 */
	void ResetBoundary(void);
	
	/*! 
	 * \overload
	 * \param[in] val_boundary - <code>TRUE</code> if the point belong to the boundary; otherwise <code>FALSE</code>.
	 */		
	void SetBoundary(bool val_boundary);
	
	/*! 
	 * \brief Provides information about if a point belong to the boundaries.
	 * \return <code>TRUE</code> if the point belong to the boundary; otherwise <code>FALSE</code>.
	 */		
	bool GetBoundary(void);
  
  /*!
	 * \brief Set if a point belong to the boundary.
	 * \param[in] val_boundary - <code>TRUE</code> if the point belong to the physical boundary; otherwise <code>FALSE</code>.
	 */
	void SetPhysicalBoundary(bool val_boundary);
  
  /*!
	 * \brief Set if a point belong to the boundary.
	 * \param[in] val_boundary - <code>TRUE</code> if the point belong to the physical boundary; otherwise <code>FALSE</code>.
	 */
	void SetSolidBoundary(bool val_boundary);
  
  /*!
	 * \brief Provides information about if a point belong to the physical boundaries (without MPI).
	 * \return <code>TRUE</code> if the point belong to the boundary; otherwise <code>FALSE</code>.
	 */
	bool GetPhysicalBoundary(void);
	
  /*!
	 * \brief Provides information about if a point belong to the physical boundaries (without MPI).
	 * \return <code>TRUE</code> if the point belong to the boundary; otherwise <code>FALSE</code>.
	 */
	bool GetSolidBoundary(void);
  
	/*! 
	 * \brief Set a color to the point that comes from the grid partitioning.
	 * \note Each domain has a different color.
	 * \param[in] val_color - Color of the point.
	 */
	void SetColor(unsigned short val_color);
	
	/*! 
	 * \brief Set the number of neighbor (artificial dissipation).
	 * \param[in] val_nneighbor - Number of neighbors.
	 */
	void SetnNeighbor(unsigned short val_nneighbor);

	/*! 
	 * \brief Get the number of neighbor of a point.
	 * \return Number of neighbors.
	 */
	unsigned short GetnNeighbor(void);
	
	/*! 
	 * \brief Get the color of a point, the color indicates to which subdomain the point belong to.
	 * \return Color of the point.
	 */
	unsigned short GetColor(void);
	
	/*! 
	 * \brief Get the global index in a parallel computation.
	 * \return Global index in a parallel computation.
	 */
	unsigned long GetGlobalIndex(void);
  
  /*!
	 * \brief Set the global index in a parallel computation.
	 * \return Global index in a parallel computation.
	 */
	void SetGlobalIndex(unsigned long val_globalindex);
	
	/*! 
	 * \brief Get the volume of the control volume at time n.
	 * \return Volume of the control volume at time n
	 */	
	su2double GetVolume_n(void);
	
	/*! 
	 * \brief Get the volume of the control volume at time n+1.
	 * \return Volume of the control volume at time n+1
	 */	
	su2double GetVolume_nM1(void);
	
	/*! 
	 * \brief Set the volume of the control volume at time n.
	 */
	void SetVolume_n(void);
	
	/*! 
	 * \brief Set the volume of the control volume at time n+1.
	 */
	void SetVolume_nM1(void);
	
  /*! 
	 * \brief Get the coordinates of the control volume at time n.
	 * \return Coordinates of the control volume at time n.
	 */	
	su2double* GetCoord_n(void);
	
	/*! 
	 * \brief Get the coordinates of the control volume at time n-1.
	 * \return Volume of the control volume at time n-1
	 */	
	su2double* GetCoord_n1(void);
	
  /*! 
	 * \brief Get the coordinates of the control volume at time n+1.
	 * \return Volume of the control volume at time n+1
	 */	
	su2double* GetCoord_p1(void);
  
	/*! 
	 * \brief Set the coordinates of the control volume at time n.
	 */
	void SetCoord_n(void);
	
	/*! 
	 * \brief Set the coordinates of the control volume at time n-1.
	 */
	void SetCoord_n1(void);

	/*!
	 * \brief Set the coordinates of the control volume at time n, for restart cases.
	 * \param[in] val_coord - Value of the grid coordinates at time n.
	 */
	void SetCoord_n(su2double *val_coord);

	/*!
	 * \brief Set the coordinates of the control volume at time n-1, for restart cases.
	 * \param[in] val_coord - Value of the grid coordinates at time n-1.
	 */
	void SetCoord_n1(su2double *val_coord);
  
	/*! 
	 * \brief Set the coordinates of the control volume at time n+1.
	 * \param[in] val_coord - Value of the grid coordinates at time n+1.
	 */	
	void SetCoord_p1(su2double *val_coord);
  
	/*! 
	 * \brief Set the volume of the control volume.
	 * \param[in] val_Volume - Value of the volume.
	 */
	void SetVolume(su2double val_Volume);
	
	/*! 
	 * \brief Set if a element is going to be moved on the deformation process.
	 * \param[in] val_move - true or false depending if the point will be moved.
	 */
	void SetMove(bool val_move);
	
	/*! 
	 * \brief Set the parent control volume of an agglomerated control volume.
	 * \param[in] val_parent_CV - Index of the parent control volume.
	 */	
	void SetParent_CV(unsigned long val_parent_CV);
	
	/*! 
	 * \brief Set the children control volumes of an agglomerated control volume.
	 * \param[in] val_nchildren_CV - Number of children.
	 * \param[in] val_children_CV - Index of the children control volume.
	 */	
	void SetChildren_CV(unsigned short val_nchildren_CV, unsigned long val_children_CV);
	
	/*! 
	 * \brief Get the parent control volume of an agglomerated control volume.
	 * \return Index of the parent control volume.
	 */	
	unsigned long GetParent_CV(void);
	
	/*! 
	 * \brief Get the children control volume of an agglomerated control volume.
	 * \param[in] val_nchildren_CV - Number of the children.
	 * \return Index of the parent control volume.
	 */	
	unsigned long GetChildren_CV(unsigned short val_nchildren_CV);
	
	/*! 
	 * \brief Get information about if a control volume has been agglomerated.
	 * \return <code>TRUE</code> if the point has been agglomerated; otherwise <code>FALSE</code>.
	 */	
	bool GetAgglomerate(void);
	
	/*! 
	 * \brief Get information about if the indirect neighbors can be agglomerated.
	 * \return <code>TRUE</code> if the indirect neigbors can be agglomerated; otherwise <code>FALSE</code>.
	 */	
	bool GetAgglomerate_Indirect(void);
	
	/*! 
	 * \brief Set information about if the indirect neighbors can be agglomerated.
	 * \param[in] val_agglomerate - The indirect neigbors can be agglomerated.
	 */	
	void SetAgglomerate_Indirect(bool val_agglomerate);
	
	/*! 
	 * \brief Get the number of children of an agglomerated control volume.
	 * \return Number of children control volume.
	 */	
	unsigned short GetnChildren_CV(void);
	
	/*! 
	 * \brief Set the number of children of an agglomerated control volume.
	 * \param[in] val_nchildren_CV - Number of children of the control volume.
	 */	
	void SetnChildren_CV(unsigned short val_nchildren_CV);
	
	/*! 
	 * \brief Get the value of the summed coordinates for implicit smoothing.
	 * \return Sum of coordinates at a point.
	 */	
	su2double *GetCoord_Sum(void);
	
	/*! 
	 * \brief Get the value of the old coordinates for implicit smoothing.
	 * \return Old coordinates at a point.
	 */	
	su2double *GetCoord_Old(void);
	
	/*! 
	 * \brief Get the value of the grid velocity at the point.
	 * \return Grid velocity at the point.
	 */	
	su2double *GetGridVel(void);
  
  /*!
	 * \brief Get the value of the grid velocity gradient at the point.
	 * \return Grid velocity gradient at the point.
	 */
	su2double **GetGridVel_Grad(void);
	
	/*!
	 * \brief Add the value of the coordinates to the <i>Coord_Sum</i> vector for implicit smoothing.
	 * \param[in] val_coord_sum - Value of the coordinates to add.
	 */	
	void AddCoord_Sum(su2double *val_coord_sum);
	
	/*!
	 * \brief Initialize the vector <i>Coord_Sum</i>.
	 */	
	void SetCoord_SumZero(void);
	
	/*!
	 * \brief Set the value of the vector <i>Coord_Old</i> for implicit smoothing.
	 * \param[in] val_coord_old - Value of the coordinates.
	 */	
	void SetCoord_Old(su2double *val_coord_old);
	
	/*! 
	 * \brief Set the value of the grid velocity at the point.
	 * \param[in] val_dim - Index of the coordinate.
	 * \param[in] val_gridvel - Value of the grid velocity.
	 */	
	void SetGridVel(unsigned short val_dim, su2double val_gridvel);
	
	/*! 
	 * \overload
	 * \param[in] val_gridvel - Value of the grid velocity.
	 */	
	void SetGridVel(su2double *val_gridvel);
  
  /*!
	 * \brief Set the gradient of the grid velocity.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	void SetGridVel_Grad(unsigned short val_var, unsigned short val_dim, su2double val_value);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	void SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	void SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG);

	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	void GetNormal(su2double *val_normal);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *		  definition of the function in all the derived classes).
	 */
	su2double *GetNormal(void);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	void SetNormal(su2double *val_face_normal);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	unsigned short GetnNodes(void);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	void SetZeroValues(void);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	void AddNormal(su2double *val_face_normal);

  /*!
   * \brief Set the adjoint values of the coordinates.
   * \param[in] adj_sol - The adjoint values of the coordinates.
   */
  void SetAdjointCoord(su2double *adj_coor);

  /*!
   * \brief Get the adjoint values of the coordinates.
   * \param[in] adj_sol - The adjoint values of the coordinates.
   */
  void GetAdjointCoord(su2double *adj_coor);

};

/*! 
 * \class CEdge
 * \brief Class for defining an edge.
 * \author F. Palacios
 */
class CEdge : public CDualGrid {
private:
	su2double *Coord_CG;			/*!< \brief Center-of-gravity of the element. */
	unsigned long *Nodes;		/*!< \brief Vector to store the global nodes of an element. */
	su2double *Normal;				/*!< \brief Normal al elemento y coordenadas de su centro de gravedad. */

public:
		
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_iPoint - First node of the edge.		 
	 * \param[in] val_jPoint - Second node of the edge.
	 * \param[in] val_nDim - Number of dimensions of the problem.		 
	 */
	CEdge(unsigned long val_iPoint, unsigned long val_jPoint, unsigned short val_nDim);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CEdge(void);
	
	/*! 
	 * \brief Set the center of gravity of the edge.
	 * \param[in] val_coord - Coordinates of all the nodes needed for computing the centre of gravity of an edge.
	 */
	void SetCoord_CG(su2double **val_coord);
	
	/*! 
	 * \brief Obtain the centre of gravity of the edge.
	 * \param[in] val_dim - Position to read the coordinate.		 
	 * \return Coordinate <i>val_dim</i> of the centre of gravity.
	 */
	su2double GetCG(unsigned short val_dim);
	
	/*! 
	 * \brief Get the nodes of the edge.
	 * \param[in] val_node - Position of the node that makes the edge.		 
	 * \return Index of the node that compose the edge.
	 */
	
	unsigned long GetNode(unsigned short val_node);
	
	/*! 
	 * \brief Get the number of nodes of an element.
	 * \return Number of nodes that set an edge (2).
	 */
	unsigned short GetnNodes(void);
	
	/*! 
	 * \brief Compute Volume associated to each edge.
	 * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
	 * \param[in] val_coord_FaceElem_CG - Coordinates of the centre of gravity of the face of an element.
	 * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
	 * \param[in] val_coord_Point - Coordinates of the point that form the control volume.		 
	 * \return Local volume associated to the edge.
	 */
	su2double GetVolume(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG, su2double *val_coord_Point);

	/*! 
	 * \overload
	 * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
	 * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
	 * \param[in] val_coord_Point - Coordinates of the point that form the control volume.
	 * \return Local volume associated to the edge.
	 */
	su2double GetVolume(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG, su2double *val_coord_Point);
	
	/*! 
	 * \brief Set the face that correspond to an edge.
	 * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
	 * \param[in] val_coord_FaceElem_CG - Coordinates of the centre of gravity of the face of an element.
	 * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] config - Definition of the particular problem.
	 * \return Compute the normal (dimensional) to the face that makes the control volume boundaries.
	 */
	void SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG);
	
	/*!
	 * \overload
	 * \brief Set the face that correspond to an edge.
	 * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
	 * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] config - Definition of the particular problem.
	 * \return Compute the normal (dimensional) to the face that makes the contorl volume boundaries.
	 */
	void SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG);
	
	/*! 
	 * \brief Copy the the normal vector of a face.
	 * \param[in] val_normal - Vector where the subroutine is goint to copy the normal (dimensional).
	 */
	void GetNormal(su2double *val_normal);
	
	/*! 
	 * \brief Get the normal to a face of the control volume asociated with an edge.
	 * \return Dimensional normal vector, the modulus is the area of the face.
	 */
	su2double *GetNormal(void);
	
	/*! 
	 * \brief Initialize normal vector.
	 */
	void SetZeroValues(void);

	/*! 
	 * \brief Set the normal vector.
	 * \param[in] val_face_normal - Vector to initialize the normal vector.
	 * \return Value of the normal vector.
	 */
	void SetNormal(su2double *val_face_normal);
	
	/*! 
	 * \brief Add a vector to the normal vector.
	 * \param[in] val_face_normal - Vector to add to the normal vector.
	 */
	void AddNormal(su2double *val_face_normal);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	su2double *GetCoord(void);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	void SetCoord(su2double *val_coord);

};

/*! 
 * \class CVertex
 * \brief Class for vertex definition (equivalent to edges, but for the boundaries).
 * \author F. Palacios
 */
class CVertex : public CDualGrid {
protected:
	unsigned long *Nodes;	/*!< \brief Vector to store the global nodes of an element. */
	su2double *Normal;			/*!< \brief Normal coordinates of the element and its center of gravity. */
	su2double Aux_Var;			/*!< \brief Auxiliar variable defined only on the surface. */
	su2double CartCoord[3];		/*!< \brief Vertex cartesians coordinates. */
	su2double VarCoord[3];		/*!< \brief Used for storing the coordinate variation due to a surface modification. */
	su2double *VarRot;   /*!< \brief Used for storing the rotation variation due to a surface modification. */
	long PeriodicPoint[5];			/*!< \brief Store the periodic point of a boundary (iProcessor, iPoint) */
  bool ActDisk_Perimeter;     /*!< \brief Identify nodes at the perimeter of the actuator disk */
	short Rotation_Type;			/*!< \brief Type of rotation associated with the vertex (MPI and periodic) */
	unsigned long Normal_Neighbor; /*!< \brief Index of the closest neighbor. */
	unsigned long *Donor_Points; /*!< \brief indices of donor points for interpolation across zones */
	unsigned long *Donor_Proc; /*!< \brief indices of donor processor for interpolation across zones in parallel */
  unsigned long Donor_Elem;   /*!< \brief Store the donor element for interpolation across zones/ */
  unsigned short Donor_Face;  /*!<\brief Store the donor face (w/in donor element) for interpolation across zones */
  su2double Basis_Function[3]; /*!< \brief Basis function values for interpolation across zones. */
  su2double *Donor_Coeff; /*!\brief Store a list of coefficients corresponding to the donor points. */
  unsigned short nDonor_Points; /*!\brief Number of points in Donor_Points; at least there will be one donor point (if the mesh is matching)*/
  
public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_point - Node of the vertex.
	 * \param[in] val_nDim - Number of dimensions of the problem.		
	 */
	CVertex(unsigned long val_point, unsigned short val_nDim);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CVertex(void);
	
	/*! 
	 * \brief Get the number of nodes of a vertex.
	 * \return Number of nodes that set a vertex (1).
	 */
	unsigned short GetnNodes(void);
	
	/*! 
	 * \brief Get the node of the vertex.
	 * \return Index of the node that compose the vertex.
	 */
	unsigned long GetNode(void);
	
	/*! 
	 * \brief Set the face that correspond to a vertex.
	 * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
	 * \param[in] val_coord_FaceElem_CG - Coordinates of the centre of gravity of the face of an element.
	 * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
	 * \return Compute the normal (dimensional) to the face that makes the vertex.
	 */
	void SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG);
	
	/*! 
	 * \overload
	 * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
	 * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
	 * \return Compute the normal (dimensional) to the face that makes the vertex.
	 */
	void SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG);
	
	/*! 
	 * \brief Copy the the normal vector of a face.
	 * \param[in] val_normal - Vector where the subroutine is goint to copy the normal (dimensional).
	 */
	void GetNormal(su2double *val_normal);
	
	/*! 
	 * \brief Get the normal to a face of the control volume asociated with a vertex.
	 * \return Dimensional normal vector, the modulus is the area of the face.
	 */
	su2double *GetNormal(void);
	
	/*! 
	 * \brief Initialize normal vector.
	 */
	void SetZeroValues(void);
	
	/*! 
	 * \brief Set the value of an auxiliary variable for gradient computation.
	 * \param[in] val_auxvar - Value of the auxiliar variable.
	 */
	void SetAuxVar(su2double val_auxvar);
	
	/*! 
	 * \brief Get the value of an auxiliary variable for gradient computation.
	 * \return Value of the auxiliar variable.
	 */
	su2double GetAuxVar(void);

  /*!
	 * \brief Add the value of an auxiliary variable for gradient computation.
	 * \param[in] val_auxvar - Value of the auxiliar variable.
	 */
	void AddAuxVar(su2double val_auxvar);
  
	/*! 
	 * \brief Set the normal vector.
	 * \param[in] val_face_normal - Vector to initialize the normal vector.
	 * \return Value of the normal vector.
	 */
	void SetNormal(su2double *val_face_normal);
	
	/*! 
	 * \brief Add a vector to the normal vector.
	 * \param[in] val_face_normal - Vector to add to the normal vector.
	 */
	void AddNormal(su2double *val_face_normal);
	
	/*! 
	 * \brief Set the value of the coordinate variation due to a surface modification.
	 * \param[in] val_varcoord - Variation of the coordinate.
	 */
	void SetVarCoord(su2double *val_varcoord);
	
	/*! 
	 * \brief Add the value of the coordinate variation due to a surface modification.
	 * \param[in] val_varcoord - Variation of the coordinate.
	 */
	void AddVarCoord(su2double *val_varcoord);
	
	/*! 
	 * \brief Get the value of the coordinate variation due to a surface modification.
	 * \return Variation of the coordinate.
	 */
	su2double *GetVarCoord(void);
	
	/*! 
	 * \brief Set the value of the cartesian coordinate for the vertex.
	 * \param[in] val_coord - Value of the cartesian coordinate.
	 */
	void SetCoord(su2double *val_coord);
	
	/*! 
	 * \brief Get the value of the cartesian coordinate for the vertex.
	 * \return Value of the cartesian coordinate of the vertex.
	 */
	su2double *GetCoord(void);
  
  /*!
	 * \brief Get the value of the cartesian coordinate for the vertex.
   * \param[in] val_dim - Variable of the dimension.
	 * \return Value of the cartesian coordinate of the vertex.
	 */
  su2double GetCoord(unsigned short val_dim);
	
	/*! 
	 * \brief Set the type of rotation associated to the vertex.
	 * \param[in] val_rotation_type - Value of the rotation that will be applied to the solution at the vertex
	 */
	void SetRotation_Type(short val_rotation_type);
	
	/*! 
	 * \brief Get the type of rotation associated to the vertex.
	 * \return Value of the rotation that must be applied to the solution of the vertex
	 */
	short GetRotation_Type(void);
  
	/*! 
	 * \overload
	 * \param[in] val_periodicpoint - Value of periodic point of the vertex.
	 * \param[in] val_processor - Processor where the point belong.
	 */
	void SetDonorPoint(long val_periodicpoint, long val_processor);
	
  /*!
   * \overload
   * \param[in] val_periodicpoint - Value of periodic point of the vertex.
   * \param[in] val_processor - Processor where the point belong.
   */
  void SetDonorPoint(long val_periodicpoint, long val_periodicglobalindex, long val_periodicvertex, long val_periodicmarker, long val_processor);
  
	/*! 
	 * \overload
	 * \param[in] val_periodicpoint - Value of periodic point of the vertex.
	 * \param[in] val_processor - Processor where the point belong.
	 * \param[in] val_globalindex - Global index of the donor point.
	 */
	void SetDonorPoint(long val_periodicpoint, long val_processor, long val_globalindex);
  
  /*!
   * \overload
   * \param[in] val_periodicpoint - Value of periodic point of the vertex.
   * \param[in] val_processor - Processor where the point belong.
   */
  void SetActDisk_Perimeter(bool val_actdisk_perimeter);

	/*!
	 * \brief Get the value of the periodic point of a vertex.
	 * \return Value of the periodic point of a vertex.
	 */
	long GetDonorPoint(void);
  
  /*!
   * \brief Get the value of the periodic point of a vertex.
   * \return Value of the periodic point of a vertex.
   */
  long GetDonorMarker(void);
  
  /*!
   * \brief Get the value of the periodic point of a vertex.
   * \return Value of the periodic point of a vertex.
   */
  long GetDonorVertex(void);

  /*!
   * \brief Get the value of the periodic point of a vertex.
   * \return Value of the periodic point of a vertex.
   */
  long GetDonorGlobalIndex(void);
  
  /*!
   * \brief Get the value of the periodic point of a vertex.
   * \return Value of the periodic point of a vertex.
   */
  long GetGlobalDonorPoint(void);

  /*!
	 * \brief Get the value of the periodic point of a vertex.
	 * \return Value of the periodic point of a vertex.
	 */
	long GetDonorProcessor(void);
  
	/*! 
	 * \brief Get the value of the periodic point of a vertex, and its somain
	 * \return Value of the periodic point of a vertex, and the domain.
	 */
	long *GetPeriodicPointDomain(void);	
  
  /*!
   * \brief Get the value of the periodic point of a vertex, and its somain
   * \return Value of the periodic point of a vertex, and the domain.
   */
  bool GetActDisk_Perimeter(void);

  /*!
	 * \brief Set the donor element of a vertex for interpolation across zones.
	 * \param[in] val_donorelem - donor element index.
	 */
	void SetDonorElem(long val_donorelem);
  
  /*!
	 * \brief Get the donor element of a vertex for interpolation across zones.
	 * \return Value of the donor element of a vertex.
	 */
	long GetDonorElem(void);

	/*!
   * \brief Set the donor face of a vertex for interpolation across zones.
   * \param[in] val_donorface- donor face index (w/in donor elem).
   */
  void SetDonorFace(unsigned short val_donorface);

  /*!
   * \brief Get the donor face of a vertex for interpolation across zones.
   * \return Value of the donor face index (w/in donor elem).
   */
  unsigned short GetDonorFace(void);
  
  /*!
	 * \brief Set the finite element basis functions needed for interpolation.
	 * \param[in] val_node - a node index of the owner element.
   * \param[in] val_basis - basis function value for the node.
	 */
	void SetBasisFunction(unsigned short val_node, su2double val_basis);
	
	/*!
	 * \brief Get the finite element basis functions needed for interpolation.
	 * \param[in] val_node - a node index of the owner element.
   * \return Value of the basis function for this node.
	 */
	su2double GetBasisFunction(unsigned short val_node);
	
	/*! 
	 * \brief Set the index of the closest neighbor to a point on the boundaries.
	 * \param[in] val_Normal_Neighbor - Index of the closest neighbor.
	 */
	void SetNormal_Neighbor(unsigned long val_Normal_Neighbor);
	
	/*! 
	 * \brief Get the value of the closest neighbor.
	 * \return Index of the closest neighbor.
	 */
	unsigned long GetNormal_Neighbor(void);
	
	/*!
   * \brief Increment the number of donor points by 1.
   */
	void IncrementnDonor(void);

	/*!
   * \brief Set the value of nDonor_Points
   * \param[in] nDonor - the number of donor points
   */
	void SetnDonorPoints(unsigned short nDonor);

	/*!
   * \brief Return the value of nDonor_Points
   * \return nDonor - the number of donor points
   */
	unsigned short GetnDonorPoints(void);

	/*!
   * \brief Set the coefficient value of a donor point.
   * \param[in] iDonor - Index of the donor point.
   * \param[in] val  - Value of the coefficent for point iDonor.
   */
	void SetDonorCoeff(unsigned short iDonor, su2double val);

	/*!
   * \brief Get the coefficient value of a donor point.
   * \param[in] iDonor - Index of the donor point.
   * \return  - Value of the coefficent for point iDonor.
   */
	su2double GetDonorCoeff(unsigned short iDonor);

  /*!
   * \brief Set the donor point of a vertex for interpolation across zones.
   * \param[in] val_donorpoint- donor face index (w/in donor elem).
   */
  void SetInterpDonorPoint(unsigned short val_donorindex, long val_donorpoint);

  /*!
   * \brief Get the value of the donor point of a vertex (for interpolation).
   * \return Value of the donor point of a vertex.
   */
  long GetInterpDonorPoint(unsigned short val_donorpoint);


  /*!
   * \brief Set the donor point of a vertex for interpolation across zones.
   * \param[in] val_donorpoint- donor face index (w/in donor elem).
   */
  void SetInterpDonorProcessor(unsigned short val_donorindex, long val_rank);

  /*!
   * \brief Get the value of the donor point of a vertex (for interpolation).
   * \return Value of the donor point of a vertex.
   */
  long GetInterpDonorProcessor(unsigned short val_donorindex);


	/*!
   * \brief Allocate memory based on how many donor points need to be stored.
   * Uses nDonor_Points
   */
	void	Allocate_DonorInfo(void);

	/*!
   * \brief Get the rotation variation
   * \return  - pointer to the vector defining the rotation
   */
  su2double *GetVarRot(void);

  /*!
   * \brief Set the rotation variation
   * \return  - pointer to the vector defining the rotation
   */
  void SetVarRot(su2double* val);

};

/*!
 * \class CTurboVertex
 * \brief Class for vertex definition for turbomachinery (equivalent to edges, but for the boundaries).
 * \author S. Vitale
 */
class CTurboVertex : public CVertex {
private:
  su2double *TurboNormal;			/*!< \brief Normal for computing correct turbomachinery quantities. */
  su2double Area;							/*!< \brief Value of the face area associated to the vertex */
  //	su2double PitchCoord;       /*!< \brief Value of the abscissa pitch wise */
  su2double AngularCoord;     /*!< \brief Value of the angular coordinate  */
  su2double DeltaAngularCoord;     /*!< \brief Value of the angular coordinate w.r.t. the minimum pitch point  */
  su2double RelAngularCoord; /*!< \brief Value of the angular coordinate w.r.t. the minimum pitch point  */

  unsigned long OldVertex;    /*!< \brief Value of the vertex numeration before the ordering */
  int GlobalIndex;						/*!< \brief Value of the vertex numeration after the ordering and global with respect to MPI partinioning */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_point - Node of the vertex.
   * \param[in] val_nDim - Number of dimensions of the problem.
   */
  CTurboVertex(unsigned long val_point, unsigned short val_nDim);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurboVertex(void);

  /*!
   * \brief set Normal in the turbomachinery frame of reference.
   * \param[in] val_normal - normal vector.
   */
  void SetTurboNormal(su2double *val_normal);

  /*!
   * \brief set face Area.
   * \param[in] val_area - value of the face area.
   */
  void SetArea(su2double val_area);

  /*!
   * \brief get face Area associate to the vertex.
   */
  su2double GetArea(void);

  /*!
   * \brief Copy the the turbo normal vector of a face.
   * \param[in] val_normal - Vector where the subroutine is goint to copy the normal (dimensionaless).
   */
  void GetTurboNormal(su2double *val_normal);

  /*!
   * \brief Get the turbo normal to a face where turboperformance are computed .
   * \return Dimensionaless normal vector, the modulus is the area of the face.
   */
  su2double *GetTurboNormal(void);

  /*!
   * \brief set vertex value not ordered.
   * \param[in] val_vertex - value of the vertex before ordering.
   */
  void SetOldVertex(unsigned long val_vertex);

  /*!
   * \brief retrieve vertex value not ordered.
   */
  unsigned long GetOldVertex(void);

  /*!
   * \brief set global index for ordered span-wise turbovertex.
   */
  void SetGlobalVertexIndex(int globalindex);

  /*!
   * \brief get global index for ordered span-wise turbovertex.
   */
  int GetGlobalVertexIndex(void);

  /*!
   * \brief set angular coord.
   */
  void SetAngularCoord(su2double angCoord);

  /*!
   * \brief get angular coord.
   */
  su2double GetAngularCoord(void);

  /*!
   * \brief set angular coord.
   */
  void SetDeltaAngularCoord(su2double deltaAngCoord);

  /*!
   * \brief get angular coord.
   */
  su2double GetDeltaAngularCoord(void);

  /*!
   * \brief set angular coord.
   */
  void SetRelAngularCoord(su2double minAngCoord);

  /*!
   * \brief get angular coord.
   */
  su2double GetRelAngularCoord(void);

};

#include "dual_grid_structure.inl"
