/*!
 * \file dual_grid_structure.hpp
 * \brief Headers of the main subroutines for doing the complete dual grid structure.
 *        The subroutines and functions are in the <i>dual_grid_structure.cpp</i> file.
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

#include <cmath>
#include <iostream>
#include <vector>

#include "config_structure.hpp"

using namespace std;

/*! 
 * \class CDualGrid
 * \brief Class for controlling the dual volume definition. The dual volume is compose by 
 *        three main elements: points, edges, and vertices.
 * \author F. Palacios.
 * \version 2.0.6
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
	~CDualGrid(void);
	
	/*! 
	 * \brief A pure virtual member.
	 */
	virtual double *GetCoord(void) = 0;
	
	/*! 
	 * \brief A pure virtual member.
	 * \param[in] val_coord - Coordinate of the point.		 
	 */
	virtual void SetCoord(double *val_coord) = 0;
	
	/*! 
	 * \brief A pure virtual member.
	 * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
	 * \param[in] val_coord_FaceElem_CG - Coordinates of the centre of gravity of the face of an element.
	 * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetNodes_Coord(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG,double *val_coord_Elem_CG, CConfig *config) = 0;
	
	/*! 
	 * \overload
	 * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
	 * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetNodes_Coord(double *val_coord_Edge_CG, double *val_coord_Elem_CG, CConfig *config) = 0;
	
	/*! 
	 * \brief A pure virtual member.
	 * \param[in] val_normal - Coordinates of the normal.
	 */
	virtual void GetNormal(double *val_normal) = 0;
	
	/*! 
	 * \brief A pure virtual member.
	 */
	virtual double *GetNormal(void) = 0;
	
	/*! 
	 * \brief A pure virtual member.
	 * \param[in] val_face_normal - Coordinates of the normal.
	 */	
	virtual void SetNormal(double *val_face_normal) = 0;

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
	virtual void AddNormal(double *val_face_normal) = 0;
};

/*! 
 * \class CPoint
 * \brief Class for point definition (including control volume definition).
 * \author F. Palacios.
 * \version 2.0.6
 */
class CPoint : public CDualGrid {
private:
	unsigned short nElem,	/*!< \brief Number of elements that set up the control volume. */
	nPoint;					/*!< \brief Number of points that set up the control volume  */
	vector<long> Elem;		/*!< \brief Elements that set up a control volume around a node. */
	vector<unsigned long> Point;	/*!< \brief Points surrounding the central node of the control volume. */
	vector<long> Edge;		/*!< \brief Edges that set up a control volume. */
	double *Volume;	/*!< \brief Volume or Area of the control volume in 3D and 2D. */
	bool Domain,		/*!< \brief Indicates if a point must be computed or belong to another boundary */
	Boundary,       /*!< \brief To see if a point belong to the boundary (including MPI). */
  PhysicalBoundary;			/*!< \brief To see if a point belong to the physical boundary (without includin MPI). */
	long *vertex; /*!< \brief Index of the vertex that correspond which the control volume (we need one for each marker in the same node). */
	double *coord,	/*!< \brief vector with the coordinates of the node. */
	*Coord_old,		/*!< \brief Old coordinates vector for geometry smoothing. */
	*Coord_sum,		/*!< \brief Sum of coordinates vector for geometry smoothing. */
  *Coord_n,		/*!< \brief Coordinates at time n for use with dynamic meshes. */
  *Coord_n1,		/*!< \brief Coordinates at time n-1 for use with dynamic meshes. */
  *Coord_p1;		/*!< \brief Coordinates at time n+1 for use with dynamic meshes. */
	double *gridvel;	/*!< \brief Velocity of the grid, in case the grid is moving. */
  double **gridvel_grad;  /*!< \brief Gradient of the grid velocity for dynamic meshes. */
	double *rotvel;         /*!< \brief Rotational velocity for a rotating frame. */
  double **rotvel_grad;   /*!< \brief Gradient of the grid velocity for a rotating frame. */
	unsigned long Parent_CV;			/*!< \brief Index of the parent control volume in the agglomeration process. */
	unsigned short nChildren_CV;		/*!< \brief Number of children in the agglomeration process. */
	vector<unsigned long> Children_CV;		/*!< \brief Index of the children control volumes in the agglomeration process. */
	bool Agglomerate_Indirect;					/*!< \brief This flag indicates if the indirect points can be agglomerated. */
	bool Agglomerate;					/*!< \brief This flag indicates if the element has been agglomerated. */
	bool Move;					/*!< \brief This flag indicates if the point is going to be move in the grid deformation process. */
	unsigned short color;	/*!< \brief Color of the point in the partitioning strategy. */
	double WallDistance;	/*!< \brief Distance to the nearest wall. */
	unsigned long GlobalIndex;	/*!< \brief Global index in the parallel simulation. */
	unsigned short nNeighbor;	/*!< \brief Color of the point in the partitioning strategy. */

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
	CPoint(double val_coord_0, double val_coord_1, unsigned long val_globalindex, CConfig *config);
	
	/*! 
	 * \overload
	 * \param[in] val_coord_0 First coordinate of the point.		 
	 * \param[in] val_coord_1 Second coordinate of the point.		
	 * \param[in] val_coord_2 Third coordinate of the point.
	 * \param[in] val_globalindex Global index in the parallel simulation.
	 * \param[in] config - Definition of the particular problem.
	 */
	CPoint(double val_coord_0, double val_coord_1, double val_coord_2, unsigned long val_globalindex, CConfig *config);
	
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
	void SetWallDistance(double val_distance);
	
	/*! 
	 * \brief Get the value of the distance to the nearest wall.
	 * \return Value of the distance to the nearest wall.
	 */
	double GetWallDistance(void);
	
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
	double GetCoord(unsigned short val_dim);
	
	/*! 
	 * \brief Get the coordinates of the control volume.
	 * \return pointer to the coordinate of the point.
	 */
	double *GetCoord(void);
	
	/*! 
	 * \brief Set the coordinates for the control volume.
	 * \param[in] val_dim - Position to store the coordinate.		 
	 * \param[in] val_coord - Coordinate for val_dim.			 
	 */
	
	void SetCoord(unsigned short val_dim, double val_coord);
	/*! 
	 * \overload
	 * \param[in] val_coord - Coordinate of the point.		 
	 */
	void SetCoord(double *val_coord);
	
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
	void AddVolume(double val_Volume);
	
	/*! 
	 * \brief Get area or volume of the control volume.
	 * \return Area or volume of the control volume.
	 */
	double GetVolume(void);
	
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
	 * \brief Provides information about if a point belong to the physical boundaries (without MPI).
	 * \return <code>TRUE</code> if the point belong to the boundary; otherwise <code>FALSE</code>.
	 */
	bool GetPhysicalBoundary(void);
	
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
	 * \brief Get the volume of the control volume at time n.
	 * \return Volume of the control volume at time n
	 */	
	double GetVolume_n(void);
	
	/*! 
	 * \brief Get the volume of the control volume at time n+1.
	 * \return Volume of the control volume at time n+1
	 */	
	double GetVolume_nM1(void);
	
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
	double* GetCoord_n(void);
	
	/*! 
	 * \brief Get the coordinates of the control volume at time n-1.
	 * \return Volume of the control volume at time n-1
	 */	
	double* GetCoord_n1(void);
	
  /*! 
	 * \brief Get the coordinates of the control volume at time n+1.
	 * \return Volume of the control volume at time n+1
	 */	
	double* GetCoord_p1(void);
  
	/*! 
	 * \brief Set the coordinates of the control volume at time n.
	 */
	void SetCoord_n(void);
	
	/*! 
	 * \brief Set the coordinates of the control volume at time n-1.
	 */
	void SetCoord_n1(void);
  
	/*! 
	 * \brief Set the coordinates of the control volume at time n+1.
	 * \param[in] val_coord - Value of the grid coordinates at time n+1.
	 */	
	void SetCoord_p1(double *val_coord);
  
	/*! 
	 * \brief Set the volume of the control volume.
	 * \param[in] val_Volume - Value of the volume.
	 */
	void SetVolume(double val_Volume);
	
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
	double *GetCoord_Sum(void);
	
	/*! 
	 * \brief Get the value of the old coordinates for implicit smoothing.
	 * \return Old coordinates at a point.
	 */	
	double *GetCoord_Old(void);
	
	/*! 
	 * \brief Get the value of the grid velocity at the point.
	 * \return Grid velocity at the point.
	 */	
	double *GetGridVel(void);
  
  /*!
	 * \brief Get the value of the grid velocity gradient at the point.
	 * \return Grid velocity gradient at the point.
	 */
	double **GetGridVel_Grad(void);
	
	/*! 
	 * \brief Get the value of the rotational velocity at the point.
	 * \return Rotational velocity at the point.
	 */	
	double *GetRotVel(void);
  
  /*!
	 * \brief Get the value of the rotational velocity gradient at the point.
	 * \return Rotational velocity gradient at the point.
	 */
	double **GetRotVel_Grad(void);
	
	/*! 
	 * \brief Add the value of the coordinates to the <i>Coord_sum</i> vector for implicit smoothing.
	 * \param[in] val_coord_sum - Value of the coordinates to add.
	 */	
	void AddCoord_Sum(double *val_coord_sum);
	
	/*! 
	 * \brief Initialize the vector <i>Coord_sum</i>.
	 */	
	void SetCoord_SumZero(void);
	
	/*! 
	 * \brief Set the value of the vector <i>Coord_old</i> for implicit smoothing.
	 * \param[in] val_coord_old - Value of the coordinates.
	 */	
	void SetCoord_Old(double *val_coord_old);
	
	/*! 
	 * \brief Set the value of the grid velocity at the point.
	 * \param[in] val_dim - Index of the coordinate.
	 * \param[in] val_gridvel - Value of the grid velocity.
	 */	
	void SetGridVel(unsigned short val_dim, double val_gridvel);
	
	/*! 
	 * \overload
	 * \param[in] val_gridvel - Value of the grid velocity.
	 */	
	void SetGridVel(double *val_gridvel);
  
  /*!
	 * \brief Set the gradient of the grid velocity.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	void SetGridVel_Grad(unsigned short val_var, unsigned short val_dim, double val_value);
	
	/*! 
	 * \brief Set the value of the rotational velocity at the point.
	 * \param[in] val_dim - Index of the coordinate.
	 * \param[in] val_rotvel - Value of the rotational velocity.
	 */	
	void SetRotVel(unsigned short val_dim, double val_rotvel);
	
	/*! 
	 * \overload
	 * \param[in] val_rotvel - Value of the rotational velocity.
	 */	
	void SetRotVel(double *val_rotvel);
  
  /*!
	 * \brief Set the gradient of the rotational grid velocity.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	void SetRotVel_Grad(unsigned short val_var, unsigned short val_dim, double val_value);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	void SetNodes_Coord(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG, double *val_coord_Elem_CG, CConfig *config);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	void SetNodes_Coord(double *val_coord_Edge_CG, double *val_coord_Elem_CG, CConfig *config);

	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	void GetNormal(double *val_normal);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *		  definition of the function in all the derived classes).
	 */
	double *GetNormal(void);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	void SetNormal(double *val_face_normal);
	
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
	void AddNormal(double *val_face_normal);
};

/*! 
 * \class CEdge
 * \brief Class for defining an edge.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CEdge : public CDualGrid {
private:
	double *Coord_CG;			/*!< \brief Center-of-gravity of the element. */
	unsigned long *Nodes;		/*!< \brief Vector to store the global nodes of an element. */
	double *Normal;				/*!< \brief Normal al elemento y coordenadas de su centro de gravedad. */
  double Rot_Flux;     /*!< \brief The exactly integrated rotational volume flux. */

public:
		
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_iPoint - First node of the edge.		 
	 * \param[in] val_jPoint - Second node of the edge.
	 * \param[in] val_ndim - Number of dimensions of the problem.		 
	 */
	CEdge(unsigned long val_iPoint, unsigned long val_jPoint, unsigned short val_ndim);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CEdge(void);
	
	/*! 
	 * \brief Set the center of gravity of the edge.
	 * \param[in] val_coord - Coordinates of all the nodes needed for computing the centre of gravity of an edge.
	 */
	void SetCG(double **val_coord);
	
	/*! 
	 * \brief Obtain the centre of gravity of the edge.
	 * \param[in] val_dim - Position to read the coordinate.		 
	 * \return Coordinate <i>val_dim</i> of the centre of gravity.
	 */
	double GetCG(unsigned short val_dim);
	
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
	double GetVolume(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG, double *val_coord_Elem_CG, double *val_coord_Point);

	/*! 
	 * \overload
	 * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
	 * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
	 * \param[in] val_coord_Point - Coordinates of the point that form the control volume.
	 * \return Local volume associated to the edge.
	 */
	double GetVolume(double *val_coord_Edge_CG, double *val_coord_Elem_CG, double *val_coord_Point);
	
	/*! 
	 * \brief Set the face that correspond to an edge.
	 * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
	 * \param[in] val_coord_FaceElem_CG - Coordinates of the centre of gravity of the face of an element.
	 * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] config - Definition of the particular problem.
	 * \return Compute the normal (dimensional) to the face that makes the control volume boundaries.
	 */
	void SetNodes_Coord(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG, double *val_coord_Elem_CG, CConfig *config);
	
	/*!
	 * \overload
	 * \brief Set the face that correspond to an edge.
	 * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
	 * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] config - Definition of the particular problem.
	 * \return Compute the normal (dimensional) to the face that makes the contorl volume boundaries.
	 */
	void SetNodes_Coord(double *val_coord_Edge_CG, double *val_coord_Elem_CG, CConfig *config);
	
	/*! 
	 * \brief Copy the the normal vector of a face.
	 * \param[in] val_normal - Vector where the subroutine is goint to copy the normal (dimensional).
	 */
	void GetNormal(double *val_normal);
	
	/*! 
	 * \brief Get the normal to a face of the control volume asociated with an edge.
	 * \return Dimensional normal vector, the modulus is the area of the face.
	 */
	double *GetNormal(void);
	
	/*! 
	 * \brief Initialize normal vector.
	 */
	void SetZeroValues(void);

	/*! 
	 * \brief Set the normal vector.
	 * \param[in] val_face_normal - Vector to initialize the normal vector.
	 * \return Value of the normal vector.
	 */
	void SetNormal(double *val_face_normal);
	
	/*! 
	 * \brief Add a vector to the normal vector.
	 * \param[in] val_face_normal - Vector to add to the normal vector.
	 */
	void AddNormal(double *val_face_normal);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	double *GetCoord(void);
	
	/*! 
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the 
	 *        definition of the function in all the derived classes).
	 */
	void SetCoord(double *val_coord);
  
  /*! 
	 * \brief Get the exact integral of the rotational volume flux.
	 * \return Value of the exactly integrated rotational volume flux.
	 */
	double GetRotFlux(void);
  
  /*! 
	 * \brief Add contribution to the exact integral of the rotational volume flux.
	 */
	void AddRotFlux(double val_rot_flux);
};

/*! 
 * \class CVertex
 * \brief Class for vertex definition (equivalent to edges, but for the boundaries).
 * \author F. Palacios.
 * \version 2.0.6
 */
class CVertex : public CDualGrid {
private:
	unsigned long *Nodes;	/*!< \brief Vector to store the global nodes of an element. */
	double *Normal;			/*!< \brief Normal al elemento y coordenadas de su centro de gravedad. */
	double Aux_Var;			/*!< \brief Auxiliar variable defined only on the surface. */
	double CarCoord[3];		/*!< \brief Vertex cartesians coordinates. */
	double VarCoord[3];		/*!< \brief Used for storing the coordinate variation due to a surface modification. */
	long PeriodicPoint[2];			/*!< \brief Store the periodic point of a boundary (iProcessor, iPoint) */
	short Rotation_Type;			/*!< \brief Type of rotation associated with the vertex (MPI and periodic) */
  short Matching_Zone;			/*!< \brief Donor zone associated with the vertex (MPI and sliding) */
  double Rot_Flux;     /*!< \brief The exactly integrated rotational volume flux. */
  bool Sharp_Corner;     /*!< \brief Flag to mark vertices at sharp corners of the surfaces. */
	unsigned long Normal_Neighbor; /*!< \brief Index of the closest neighbor. */
  unsigned long Donor_Elem;   /*!< \brief Store the donor element for interpolation across zones/ */
  double Basis_Function[3]; /*!< \brief Basis function values for interpolation across zones. */
  
public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_point - Node of the vertex.
	 * \param[in] val_ndim - Number of dimensions of the problem.		
	 */
	CVertex(unsigned long val_point, unsigned short val_ndim);
	
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
   * \param[in] config - Definition of the particular problem.
	 * \return Compute the normal (dimensional) to the face that makes the vertex.
	 */
	void SetNodes_Coord(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG, double *val_coord_Elem_CG, CConfig *config);
	
	/*! 
	 * \overload
	 * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
	 * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] config - Definition of the particular problem.
	 * \return Compute the normal (dimensional) to the face that makes the vertex.
	 */
	void SetNodes_Coord(double *val_coord_Edge_CG, double *val_coord_Elem_CG, CConfig *config);
	
	/*! 
	 * \brief Copy the the normal vector of a face.
	 * \param[in] val_normal - Vector where the subroutine is goint to copy the normal (dimensional).
	 */
	void GetNormal(double *val_normal);
	
	/*! 
	 * \brief Get the normal to a face of the control volume asociated with a vertex.
	 * \return Dimensional normal vector, the modulus is the area of the face.
	 */
	double *GetNormal(void);
	
	/*! 
	 * \brief Initialize normal vector.
	 */
	void SetZeroValues(void);
	
	/*! 
	 * \brief Set the value of an auxiliary variable for gradient computation.
	 * \param[in] val_auxvar - Value of the auxiliar variable.
	 */
	void SetAuxVar(double val_auxvar);
	
	/*! 
	 * \brief Get the value of an auxiliary variable for gradient computation.
	 * \return Value of the auxiliar variable.
	 */
	double GetAuxVar(void);

  /*!
	 * \brief Add the value of an auxiliary variable for gradient computation.
	 * \param[in] val_auxvar - Value of the auxiliar variable.
	 */
	void AddAuxVar(double val_auxvar);
  
	/*! 
	 * \brief Set the normal vector.
	 * \param[in] val_face_normal - Vector to initialize the normal vector.
	 * \return Value of the normal vector.
	 */
	void SetNormal(double *val_face_normal);
	
	/*! 
	 * \brief Add a vector to the normal vector.
	 * \param[in] val_face_normal - Vector to add to the normal vector.
	 */
	void AddNormal(double *val_face_normal);
	
	/*! 
	 * \brief Set the value of the coordinate variation due to a surface modification.
	 * \param[in] val_varcoord - Variation of the coordinate.
	 */
	void SetVarCoord(double *val_varcoord);
	
	/*! 
	 * \brief Add the value of the coordinate variation due to a surface modification.
	 * \param[in] val_varcoord - Variation of the coordinate.
	 */
	void AddVarCoord(double *val_varcoord);
	
	/*! 
	 * \brief Get the value of the coordinate variation due to a surface modification.
	 * \return Variation of the coordinate.
	 */
	double *GetVarCoord(void);
	
	/*! 
	 * \brief Set the value of the cartesian coordinate for the vertex.
	 * \param[in] val_coord - Value of the cartesian coordinate.
	 */
	void SetCoord(double *val_coord);
	
	/*! 
	 * \brief Get the value of the cartesian coordinate for the vertex.
	 * \return Value of the cartesian coordinate of the vertex.
	 */
	double *GetCoord(void);
	
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
	 * \brief Get the matching zone index for a sliding interface vertex.
	 * \return Matching zone for a sliding interface vertex.
	 */
  short GetMatching_Zone(void);
	
	/*!
	 * \brief Set the matching zone index for a sliding interface vertex.
	 * \param[in] val_matching_zone - Matching zone index for a sliding interface vertex.
	 */
  void SetMatching_Zone(short val_matching_zone);
	
	/*! 
	 * \brief Set the periodic point of a vertex.
	 * \param[in] val_periodicpoint - Value of periodic point of the vertex.
	 */
	void SetDonorPoint(long val_periodicpoint);
  
	/*! 
	 * \overload
	 * \param[in] val_periodicpoint - Value of periodic point of the vertex.
	 * \param[in] val_processor - Processor where the point belong.
	 */
	void SetDonorPoint(long val_periodicpoint, long val_processor);
	
	/*! 
	 * \brief Get the value of the periodic point of a vertex.
	 * \return Value of the periodic point of a vertex.
	 */
	long GetDonorPoint(void);	
  
	/*! 
	 * \brief Get the value of the periodic point of a vertex, and its somain
	 * \return Value of the periodic point of a vertex, and the domain.
	 */
	long *GetPeriodicPointDomain(void);	
  
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
	 * \brief Set the finite element basis functions needed for interpolation.
	 * \param[in] val_node - a node index of the owner element.
   * \param[in] val_basis - basis function value for the node.
	 */
	void SetBasisFunction(unsigned short val_node, double val_basis);
	
	/*!
	 * \brief Get the finite element basis functions needed for interpolation.
	 * \param[in] val_node - a node index of the owner element.
   * \return Value of the basis function for this node.
	 */
	double GetBasisFunction(unsigned short val_node);
  
  /*! 
	 * \brief Get the exact integral of the rotational volume flux.
	 * \return Value of the exactly integrated rotational volume flux.
	 */
	double GetRotFlux(void);
  
  /*! 
	 * \brief Add contribution to the exact integral of the rotational volume flux.
	 */
	void AddRotFlux(double val_rot_flux);
  
  /*! 
	 * \brief Set the boolean for a corner vertex.
	 * \param[in] val_sharp_corner - <code>TRUE</code> if this vertex sits on a sharp corner; otherwise <code>FALSE</code>.
	 */
	void SetSharp_Corner(bool val_sharp_corner);
	
	/*! 
	 * \brief Get the value of an auxiliar variable for gradient computation.
	 * \return <code>TRUE</code> if this vertex sits on a sharp corner; otherwise <code>FALSE</code>.
	 */
	bool GetSharp_Corner(void);
	
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
	
};

#include "dual_grid_structure.inl"
