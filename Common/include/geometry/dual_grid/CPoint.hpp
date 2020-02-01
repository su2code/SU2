/*!
 * \file CPoint.hpp
 * \brief Headers of the main subroutines for doing the complete dual grid structure.
 *        The subroutines and functions are in the <i>CPoint.cpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 7.0.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "CDualGrid.hpp"

/*!
 * \class CPoint
 * \brief Class for point definition (including control volume definition).
 * \author F. Palacios
 */
class CPoint final : public CDualGrid {
private:
  unsigned short nElem,               /*!< \brief Number of elements that set up the control volume. */
  nPoint;                             /*!< \brief Number of points that set up the control volume  */
  vector<long> Elem;                  /*!< \brief Elements that set up a control volume around a node. */
  vector<unsigned long> Point;        /*!< \brief Points surrounding the central node of the control volume. */
  vector<long> Edge;                  /*!< \brief Edges that set up a control volume. */
  su2double *Volume;                  /*!< \brief Volume or Area of the control volume in 3D and 2D. */
  su2double Periodic_Volume;          /*!< \brief Missing component of volume or area of a control volume on a periodic marker in 3D and 2D. */
  bool Domain,                        /*!< \brief Indicates if a point must be computed or belong to another boundary */
  Boundary,                           /*!< \brief To see if a point belong to the boundary (including MPI). */
  PhysicalBoundary,                   /*!< \brief To see if a point belong to the physical boundary (without includin MPI). */
  SolidBoundary,                      /*!< \brief To see if a point belong to the physical boundary (without includin MPI). */
  PeriodicBoundary;                   /*!< \brief To see if a point belongs to a periodic boundary (without including MPI). */
  long *Vertex;                       /*!< \brief Index of the vertex that correspond which the control volume (we need one for each marker in the same node). */
  su2double *Coord,                   /*!< \brief vector with the coordinates of the node. */
            *Coord_Old,               /*!< \brief Old coordinates vector for primal solution reloading for Disc.Adj. with dynamic grid. */
            *Coord_Sum,               /*!< \brief Sum of coordinates vector for geometry smoothing. */
            *Coord_n,                 /*!< \brief Coordinates at time n for use with dynamic meshes. */
            *Coord_n1,                /*!< \brief Coordinates at time n-1 for use with dynamic meshes. */
            *Coord_p1;                /*!< \brief Coordinates at time n+1 for use with dynamic meshes. */
  su2double *GridVel;                 /*!< \brief Velocity of the grid for dynamic mesh cases. */
  su2double **GridVel_Grad;           /*!< \brief Gradient of the grid velocity for dynamic meshes. */
  unsigned long Parent_CV;            /*!< \brief Index of the parent control volume in the agglomeration process. */
  unsigned short nChildren_CV;        /*!< \brief Number of children in the agglomeration process. */
  vector<unsigned long> Children_CV;  /*!< \brief Index of the children control volumes in the agglomeration process. */
  bool Agglomerate_Indirect,          /*!< \brief This flag indicates if the indirect points can be agglomerated. */
  Agglomerate;                        /*!< \brief This flag indicates if the element has been agglomerated. */
  bool Move;                          /*!< \brief This flag indicates if the point is going to be move in the grid deformation process. */
  unsigned long color;                /*!< \brief Color of the point in the partitioning strategy. */
  su2double Wall_Distance;            /*!< \brief Distance to the nearest wall. */
  su2double SharpEdge_Distance;       /*!< \brief Distance to a sharp edge. */
  su2double Curvature;                /*!< \brief Value of the surface curvature (SU2_GEO). */
  unsigned long GlobalIndex;          /*!< \brief Global index in the parallel simulation. */
  unsigned short nNeighbor;           /*!< \brief Number of neighbors. */
  bool Flip_Orientation;              /*!< \brief Flip the orientation of the normal. */
  su2double MaxLength;                /*!< \brief The maximum cell-center to cell-center length. */
  int *AD_InputIndex,                 /*!< \brief Indices of Coord variables in the adjoint vector. */
  *AD_OutputIndex;                    /*!< \brief Indices of Coord variables in the adjoint vector after having been updated. */

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
  ~CPoint(void) override;

  /*!
   * \brief For parallel computation, its indicates if a point must be computed or not.
   * \param[in] val_domain - <code>TRUE</code> if the point belong to the domain; otherwise <code>FALSE</code>.
   */
  inline void SetDomain(bool val_domain) { Domain = val_domain; }

  /*!
   * \brief For parallel computation, its indicates if a point must be computed or not.
   * \return <code>TRUE</code> if the node belong to the physical domain; otherwise <code>FALSE</code>.
   */
  inline bool GetDomain(void) const { return Domain; }

  /*!
   * \brief Set the value of the distance to the nearest wall.
   * \param[in] val_distance - Value of the distance.
   */
  inline void SetWall_Distance(su2double val_distance) { Wall_Distance = val_distance; }

  /*!
   * \brief Set the value of the distance to a sharp edge.
   * \param[in] val_distance - Value of the distance.
   */
  inline void SetSharpEdge_Distance(su2double val_distance) { SharpEdge_Distance = val_distance; }

  /*!
   * \brief Get the value of the distance to the nearest wall.
   * \return Value of the distance to the nearest wall.
   */
  inline su2double GetWall_Distance(void) const { return Wall_Distance; }

  /*!
   * \brief Set the value of the curvature at a surface node.
   * \param[in] val_curvature - Value of the curvature.
   */
  inline void SetCurvature(su2double val_curvature) { Curvature = val_curvature; }

  /*!
   * \brief Get the value of the curvature at a surface node.
   * \return Value of the curvature.
   */
  inline su2double GetCurvature(void) const { return Curvature; }

  /*!
   * \brief Get the value of the distance to a sharp edge
   * \return Value of the distance to the nearest wall.
   */
  inline su2double GetSharpEdge_Distance(void) const { return SharpEdge_Distance; }

  /*!
   * \brief Set the number of elements that compose the control volume.
   * \param[in] val_nElem - Number of elements that make the control volume around a node.
   */
  inline void SetnElem(unsigned short val_nElem) { nElem = val_nElem; }

  /*!
   * \brief Set the number of points that compose the control volume.
   * \param[in] val_nPoint - Number of points that compose the control volume (points surrounding points).
   */
  inline void SetnPoint(unsigned short val_nPoint) { nPoint = val_nPoint; }

  /*!
   * \brief Get the coordinates dor the control volume.
   * \param[in] val_dim - Number of dimensions of the problem.
   * \return Coordinate that correspond with <i>val_dim</i>.
   */
  inline su2double GetCoord(unsigned short val_dim) const { return Coord[val_dim]; }

  /*!
   * \brief Get the coordinates of the control volume.
   * \return pointer to the coordinate of the point.
   */
  inline su2double *GetCoord(void) override { return Coord; }

  /*!
   * \brief Set the coordinates for the control volume.
   * \param[in] val_dim - Position to store the coordinate.
   * \param[in] val_coord - Coordinate for val_dim.
   */
  inline void SetCoord(unsigned short val_dim, su2double val_coord) { Coord[val_dim] = val_coord; }

  /*!
   * \brief Set the adjoint vector indices of Coord vector.
   * \param[in] input - Save them to the input or output indices vector.
   */
  void SetIndex(bool input);

  /*!
   * \brief Set the adjoint values of the (geometric) coordinates.
   * \param[in] adj_sol - Adjoint values of the Coord variables.
   */
  void SetAdjointSolution(const su2double *adj_sol);

  /*!
   * \brief Get the adjoint values of the (geometric) coordinates.
   * \param[in] adj_sol - Adjoint values of the Coord variables.
   */
  su2double GetAdjointSolution(unsigned short iDim) const;

  /*!
   * \brief Get the coordinates of the control volume.
   * \return pointer to the coordinate of the point.
   */
  inline bool GetFlip_Orientation(void) const { return Flip_Orientation; }

  /*!
   * \brief Set the coordinates for the control volume.
   * \param[in] val_dim - Position to store the coordinate.
   * \param[in] val_coord - Coordinate for val_dim.
   */
  inline void SetFlip_Orientation(void) { Flip_Orientation = true; }

  /*!
   * \brief Set the coordinates for the control volume.
   * \param[in] val_dim - Position to store the coordinate.
   * \param[in] val_coord - Coordinate for val_dim.
   */
  inline void AddCoord(unsigned short val_dim, su2double val_coord) { Coord[val_dim] += val_coord; }

  /*!
   * \overload
   * \param[in] val_coord - Coordinate of the point.
   */
  inline void SetCoord(const su2double *val_coord) override {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Coord[iDim]=val_coord[iDim];
  }

  /*!
   * \brief Get the number of elements that compose the control volume.
   * \return Number of elements that compose the control volume.
   */
  inline unsigned short GetnElem(void) const { return nElem; }

  /*!
   * \brief Get the number of points that compose the control volume.
   * \return Number of points that compose the control volume.
   */
  inline unsigned short GetnPoint(void) const { return nPoint; }

  /*!
   * \brief Set the elements that set the control volume.
   * \param[in] val_elem - Element to be added.
   */
  inline void SetElem(unsigned long val_elem) { Elem.push_back(val_elem); nElem = Elem.size(); }

  /*!
   * \brief Reset the elements of a control volume.
   */
  inline void ResetElem(void) { Elem.clear(); nElem = 0; }

  /*!
   * \brief Reset the points that compose the control volume.
   */
  inline void ResetPoint(void) { Point.clear(); Edge.clear(); nPoint = 0; }

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
  inline void SetEdge(long val_edge, unsigned short val_nedge) { Edge[val_nedge] = val_edge; }

  /*!
   * \brief Set the boundary vertex that compose the control volume.
   * \param[in] val_vertex - Vertex to be added.
   * \param[in] val_nMarker - Marker of the vertex to be added (position where is going to be stored).
   */
  inline void SetVertex(long val_vertex, unsigned short val_nmarker) {
    if (Boundary) Vertex[val_nmarker] = val_vertex;
  }

  /*!
   * \brief Get all the elements that compose the control volume.
   * \param[in] val_elem - Position where the element is stored.
   * \return Index of the element.
   */
  inline unsigned long GetElem(unsigned short val_elem) const { return Elem[val_elem]; }

  /*!
   * \brief Get all the points that compose the control volume.
   * \param[in] val_point - Position where the point is stored.
   * \return Index of the point.
   */
  inline unsigned long GetPoint(unsigned short val_point) const { return Point[val_point]; }

  /*!
   * \brief Get all the edges that compose the control volume.
   * \param[in] val_edge - Position where the edge is stored.
   * \return Index of the edge.
   */
  inline long GetEdge(unsigned short val_edge) const { return Edge[val_edge]; }

  /*!
   * \brief Get the vertex that compose the control volume for a marker.
   * \param[in] val_marker - Position where the vertex is stored.
   * \return Index of the vertex.
   */
  inline long GetVertex(unsigned short val_marker) const {
    if (Boundary) return Vertex[val_marker];
    else return -1;
  }

  /*!
   * \brief Adds some area or volume of the CV.
   * \param[in] val_Volume - Local volume to be added to the total one.
   */
  inline void AddVolume(su2double val_Volume) { Volume[0] += val_Volume; }

  /*!
   * \brief Get area or volume of the control volume.
   * \return Area or volume of the control volume.
   */
  inline su2double GetVolume(void) const { return Volume[0]; }

  /*!
   * \brief Get the missing component of area or volume for a control volume on a periodic marker.
   * \return Periodic component of area or volume for a control volume on a periodic marker.
   */
  inline su2double GetPeriodicVolume(void) const { return Periodic_Volume; }

  /*!
   * \brief Set the missing component of area or volume for a control volume on a periodic marker.
   * \param[in] val_volume - Value of the volume from the missing components of the CV on the periodic marker.
   */
  inline void SetPeriodicVolume(su2double val_volume) { Periodic_Volume = val_volume; }

  /*!
   * \brief Get the maximum cell-center to cell-center length.
   * \return The maximum cell-center to cell-center length.
   */
  inline su2double GetMaxLength(void) const {return MaxLength;}

  /*!
   * \brief Get information about the movement of the node.
   * \return <code>TRUE</code> if the point is going to be moved; otherwise <code>FALSE</code>.
   */
  inline bool GetMove(void) const { return Move; }

  /*!
   * \brief Set if a point belong to the boundary.
   * \note It also create the structure to store the vertex.
   * \param[in] val_nmarker - Max number of marker.
   */
  void SetBoundary(unsigned short val_nmarker);

  /*!
   * \brief Reset the boundary of a control volume.
   */
  inline void ResetBoundary(void) { if (Vertex != NULL) delete [] Vertex; Boundary = false; }

  /*!
   * \overload
   * \param[in] val_boundary - <code>TRUE</code> if the point belong to the boundary; otherwise <code>FALSE</code>.
   */
  inline void SetBoundary(bool val_boundary) { Boundary = val_boundary; }

  /*!
   * \brief Provides information about if a point belong to the boundaries.
   * \return <code>TRUE</code> if the point belong to the boundary; otherwise <code>FALSE</code>.
   */
  inline bool GetBoundary(void) const { return Boundary; }

  /*!
   * \brief Set if a point belong to the boundary.
   * \param[in] val_boundary - <code>TRUE</code> if the point belong to the physical boundary; otherwise <code>FALSE</code>.
   */
  inline void SetPhysicalBoundary(bool val_boundary) { PhysicalBoundary = val_boundary; }

  /*!
   * \brief Set if a point belong to the boundary.
   * \param[in] val_boundary - <code>TRUE</code> if the point belong to the physical boundary; otherwise <code>FALSE</code>.
   */
  inline void SetSolidBoundary(bool val_boundary) { SolidBoundary = val_boundary; }

  /*!
   * \brief Set if a point belongs to a periodic boundary.
   * \param[in] val_boundary - <code>TRUE</code> if the point belongs to a periodic boundary; otherwise <code>FALSE</code>.
   */
  inline void SetPeriodicBoundary(bool val_boundary) { PeriodicBoundary = val_boundary; }

  /*!
   * \brief Provides information about if a point belong to the physical boundaries (without MPI).
   * \return <code>TRUE</code> if the point belong to the boundary; otherwise <code>FALSE</code>.
   */
  inline bool GetPhysicalBoundary(void) const { return PhysicalBoundary; }

  /*!
   * \brief Provides information about if a point belong to the physical boundaries (without MPI).
   * \return <code>TRUE</code> if the point belong to the boundary; otherwise <code>FALSE</code>.
   */
  inline bool GetSolidBoundary(void) const { return SolidBoundary; }

  /*!
   * \brief Provides information about if a point belongs to a periodic boundary (without MPI).
   * \return <code>TRUE</code> if the point belongs to a periodic boundary; otherwise <code>FALSE</code>.
   */
  inline bool GetPeriodicBoundary(void) const { return PeriodicBoundary; }

  /*!
   * \brief Set a color to the point that comes from the grid partitioning.
   * \note Each domain has a different color.
   * \param[in] val_color - Color of the point.
   */
  inline void SetColor(unsigned long val_color) { color = val_color; }

  /*!
   * \brief Set the number of neighbor (artificial dissipation).
   * \param[in] val_nneighbor - Number of neighbors.
   */
  inline void SetnNeighbor(unsigned short val_nneighbor) { nNeighbor = val_nneighbor; }

  /*!
   * \brief Get the number of neighbor of a point.
   * \return Number of neighbors.
   */
  inline unsigned short GetnNeighbor(void) const { return nNeighbor; }

  /*!
   * \brief Get the color of a point, the color indicates to which subdomain the point belong to.
   * \return Color of the point.
   */
  inline unsigned long GetColor(void) const { return color; }

  /*!
   * \brief Get the global index in a parallel computation.
   * \return Global index in a parallel computation.
   */
  inline unsigned long GetGlobalIndex(void) const { return GlobalIndex; }

  /*!
   * \brief Set the global index in a parallel computation.
   * \return Global index in a parallel computation.
   */
  inline void SetGlobalIndex(unsigned long val_globalindex) { GlobalIndex = val_globalindex; }

  /*!
   * \brief Get the volume of the control volume at time n.
   * \return Volume of the control volume at time n
   */
  inline su2double GetVolume_n(void) const { return Volume[1]; }

  /*!
   * \brief Get the volume of the control volume at time n+1.
   * \return Volume of the control volume at time n+1
   */
  inline su2double GetVolume_nM1(void) const { return Volume[2]; }

  /*!
   * \brief Set the volume of the control volume at time n.
   */
  inline void SetVolume_n(void) { Volume[1] = Volume[0]; }

  /*!
   * \brief Set the volume of the control volume at time n+1.
   */
  inline void SetVolume_nM1(void) { Volume[2] = Volume[1]; }

  /*!
   * \brief Get the coordinates of the control volume at time n.
   * \return Coordinates of the control volume at time n.
   */
  inline su2double *GetCoord_n(void) { return Coord_n; }

  /*!
   * \brief Get the coordinates of the control volume at time n-1.
   * \return Volume of the control volume at time n-1
   */
  inline su2double *GetCoord_n1(void) { return Coord_n1; }

  /*!
   * \brief Get the coordinates of the control volume at time n+1.
   * \return Volume of the control volume at time n+1
   */
  inline su2double *GetCoord_p1(void) { return Coord_p1; }

  /*!
   * \brief Set the coordinates of the control volume at time n to the ones in <i>Coord</i>.
   */
  inline void SetCoord_n(void) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Coord_n[iDim] = Coord[iDim];
  }

  /*!
   * \brief Set the coordinates of the control volume at time n-1 to the ones in <i>Coord_n</i>.
   */
  inline void SetCoord_n1(void) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Coord_n1[iDim] = Coord_n[iDim];
  }

  /*!
   * \brief Set the coordinates of the control volume at time n, for restart cases.
   * \param[in] val_coord - Value of the grid coordinates at time n.
   */
  inline void SetCoord_n(su2double *val_coord) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Coord_n[iDim] = val_coord[iDim];
  }

  /*!
   * \brief Set the coordinates of the control volume at time n-1, for restart cases.
   * \param[in] val_coord - Value of the grid coordinates at time n-1.
   */
  inline void SetCoord_n1(su2double *val_coord) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Coord_n1[iDim] = val_coord[iDim];
  }
  /*!
   * \brief Set the coordinates of the control volume at time n+1.
   * \param[in] val_coord - Value of the grid coordinates at time n+1.
   */
  inline void SetCoord_p1(su2double *val_coord) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Coord_p1[iDim] = val_coord[iDim];
  }

  /*!
   * \brief Set the volume of the control volume.
   * \param[in] val_Volume - Value of the volume.
   */
  inline void SetVolume(su2double val_Volume) { Volume[0] = val_Volume; }

  /*!
   * \brief Set the max cell-center to cell-center length.
   * \param[in] val_max_length - Value of the max length
   */
  inline void SetMaxLength(su2double val_max_length) { MaxLength = val_max_length; }

  /*!
   * \brief Set if a element is going to be moved on the deformation process.
   * \param[in] val_move - true or false depending if the point will be moved.
   */
  inline void SetMove(bool val_move) { Move = val_move; }

  /*!
   * \brief Set the parent control volume of an agglomerated control volume.
   * \param[in] val_parent_CV - Index of the parent control volume.
   */
  inline void SetParent_CV(unsigned long val_parent_CV) { Parent_CV = val_parent_CV; Agglomerate = true; }

  /*!
   * \brief Set the children control volumes of an agglomerated control volume.
   * \param[in] val_nchildren_CV - Number of children.
   * \param[in] val_children_CV - Index of the children control volume.
   */
  inline void SetChildren_CV(unsigned short val_nchildren_CV, unsigned long val_children_CV) {
    if (Children_CV.size() <= val_nchildren_CV) Children_CV.resize(val_nchildren_CV+1);
    Children_CV[val_nchildren_CV] = val_children_CV;
  }

  /*!
   * \brief Get the parent control volume of an agglomerated control volume.
   * \return Index of the parent control volume.
   */
  inline unsigned long GetParent_CV(void) const { return Parent_CV; }

  /*!
   * \brief Get the children control volume of an agglomerated control volume.
   * \param[in] val_nchildren_CV - Number of the children.
   * \return Index of the parent control volume.
   */
  inline unsigned long GetChildren_CV(unsigned short val_nchildren_CV) const {  return Children_CV[val_nchildren_CV]; }

  /*!
   * \brief Get information about if a control volume has been agglomerated.
   * \return <code>TRUE</code> if the point has been agglomerated; otherwise <code>FALSE</code>.
   */
  inline bool GetAgglomerate(void) const { return Agglomerate; }

  /*!
   * \brief Get information about if the indirect neighbors can be agglomerated.
   * \return <code>TRUE</code> if the indirect neigbors can be agglomerated; otherwise <code>FALSE</code>.
   */
  inline bool GetAgglomerate_Indirect(void) const { return Agglomerate_Indirect; }

  /*!
   * \brief Set information about if the indirect neighbors can be agglomerated.
   * \param[in] val_agglomerate - The indirect neigbors can be agglomerated.
   */
  inline void SetAgglomerate_Indirect(bool val_agglomerate) { Agglomerate_Indirect = val_agglomerate; };

  /*!
   * \brief Get the number of children of an agglomerated control volume.
   * \return Number of children control volume.
   */
  inline unsigned short GetnChildren_CV(void) const { return nChildren_CV; }

  /*!
   * \brief Set the number of children of an agglomerated control volume.
   * \param[in] val_nchildren_CV - Number of children of the control volume.
   */
  inline void SetnChildren_CV(unsigned short val_nchildren_CV) {  nChildren_CV = val_nchildren_CV; }

  /*!
   * \brief Get the value of the summed coordinates for implicit smoothing.
   * \return Sum of coordinates at a point.
   */
  inline su2double *GetCoord_Sum(void) { return Coord_Sum; }

  /*!
   * \brief Get the value of the old coordinates for implicit smoothing.
   * \return Old coordinates at a point.
   */
  inline su2double *GetCoord_Old(void) { return Coord_Old; }

  /*!
   * \brief Get the value of the grid velocity at the point.
   * \return Grid velocity at the point.
   */
  inline su2double *GetGridVel(void) { return GridVel; }

  /*!
   * \brief Get the value of the grid velocity gradient at the point.
   * \return Grid velocity gradient at the point.
   */
  inline su2double **GetGridVel_Grad(void) { return GridVel_Grad; }

  /*!
   * \brief Add the value of the coordinates to the <i>Coord_Sum</i> vector for implicit smoothing.
   * \param[in] val_coord_sum - Value of the coordinates to add.
   */
  inline void AddCoord_Sum(const su2double *val_coord_sum) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Coord_Sum[iDim] += val_coord_sum[iDim];
  }

  /*!
   * \brief Initialize the vector <i>Coord_Sum</i>.
   */
  inline void SetCoord_SumZero(void) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Coord_Sum[iDim] = 0.0;
  }

  /*!
   * \brief Set the value of the vector <i>Coord_Old</i> for implicit smoothing.
   * \param[in] val_coord_old - Value of the coordinates.
   */
  inline void SetCoord_Old(const su2double *val_coord_old) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Coord_Old[iDim] = val_coord_old[iDim];
  }

  /*!
   * \brief Set the value of the vector <i>Coord_Old</i> to <i>Coord</i>.
   */
  inline void SetCoord_Old (void) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Coord_Old[iDim] = Coord[iDim];
  }

  /*!
   * \brief Set the value of the grid velocity at the point.
   * \param[in] val_dim - Index of the coordinate.
   * \param[in] val_gridvel - Value of the grid velocity.
   */
  inline void SetGridVel(unsigned short val_dim, su2double val_gridvel) { GridVel[val_dim] = val_gridvel; }

  /*!
   * \overload
   * \brief Set the value of the grid velocity at the point.
   * \param[in] val_gridvel - Value of the grid velocity.
   */
  inline void SetGridVel(const su2double *val_gridvel) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      GridVel[iDim] = val_gridvel[iDim];
  }

  /*!
   * \brief Set the gradient of the grid velocity.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value of the gradient.
   */
  inline void SetGridVel_Grad(unsigned short val_var, unsigned short val_dim, su2double val_value) {
    GridVel_Grad[val_var][val_dim] = val_value;
  }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline void SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG,
                             su2double *val_coord_Elem_CG) override { }


  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline void SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG) override { }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline void GetNormal(su2double *val_normal) const override { }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *      definition of the function in all the derived classes).
   */
  inline su2double *GetNormal(void) override { return nullptr; }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline void SetNormal(const su2double *val_face_normal) override { }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline unsigned short GetnNodes() const override { return 0; }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline void SetZeroValues(void) override { }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline void AddNormal(const su2double *val_face_normal) override { }

  /*!
   * \brief Set the adjoint values of the coordinates.
   * \param[in] adj_sol - The adjoint values of the coordinates.
   */
  inline void SetAdjointCoord(const su2double *adj_coor){
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      SU2_TYPE::SetDerivative(Coord[iDim], SU2_TYPE::GetValue(adj_coor[iDim]));
  }

  /*!
   * \brief Set the adjoint values of the coordinates.
   * \param[in] adj_sol - The adjoint values of the coordinates.
   */
  inline void SetAdjointCoord_LocalIndex(const su2double *adj_coor){
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      AD::SetDerivative(AD_OutputIndex[iDim], SU2_TYPE::GetValue(adj_coor[iDim]));
  }

  /*!
   * \brief Get the adjoint values of the coordinates.
   * \param[in] adj_sol - The adjoint values of the coordinates.
   */
  inline void GetAdjointCoord(su2double *adj_coor) const {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      adj_coor[iDim] = SU2_TYPE::GetDerivative(Coord[iDim]);
  }

  /*!
   * \brief Get the adjoint values of the coordinates.
   * \param[in] adj_sol - The adjoint values of the coordinates.
   */
  inline void GetAdjointCoord_LocalIndex(su2double *adj_coor) const {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      adj_coor[iDim] = AD::GetDerivative(AD_InputIndex[iDim]);
  }

};
