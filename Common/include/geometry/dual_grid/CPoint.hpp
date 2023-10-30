/*!
 * \file CPoint.hpp
 * \brief Declaration of the point class that stores geometric and adjacency
 *        information for dual control volumes.
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../containers/C2DContainer.hpp"
#include "../../containers/container_decorators.hpp"
#include "../../toolboxes/graph_toolbox.hpp"
#include <vector>
#include "../../toolboxes/ndflattener.hpp"

using namespace std;

class CConfig;
class CPhysicalGeometry;

/*!
 * \class CPoint
 * \brief Class for point definition (dual control volumes).
 * \author F. Palacios
 */
class CPoint {
 private:
  friend class CPhysicalGeometry;

  const unsigned long nDim = 0;

  su2vector<unsigned long> GlobalIndex; /*!< \brief Global index in the parallel simulation. */
  su2vector<unsigned long> Color;       /*!< \brief Color of the point in the partitioning strategy. */

  CCompressedSparsePatternUL Point; /*!< \brief Points surrounding the central node of the control volume. */
  CCompressedSparsePatternL Edge;   /*!< \brief Edges that set up a control volume (same sparse structure as Point). */
  CCompressedSparsePatternL Elem;   /*!< \brief Elements that set up a control volume around a node. */
  vector<vector<long> > Vertex; /*!< \brief Index of the vertex that correspond which the control volume (we need one
                                   for each marker in the same node). */

  su2activevector Volume;          /*!< \brief Volume or Area of the control volume in 3D and 2D. */
  su2activevector Volume_n;        /*!< \brief Volume at time n. */
  su2activevector Volume_nM1;      /*!< \brief Volume at time n-1. */
  su2activevector Volume_Old;      /*!< \brief Old containers for Volume. */
  su2activevector Volume_n_Old;    /*!< \brief Old containers for Volume at time n. */
  su2activevector Volume_nM1_Old;  /*!< \brief Old containers for Volume at time n-1. */
  su2activevector Periodic_Volume; /*!< \brief Missing component of volume or area of a control volume on a periodic
                                      marker in 3D and 2D. */

  su2vector<bool> Domain;   /*!< \brief Indicates if a point must be computed or belong to another boundary */
  su2vector<bool> Boundary; /*!< \brief To see if a point belong to the boundary (including MPI). */
  su2vector<bool>
      PhysicalBoundary; /*!< \brief To see if a point belong to the physical boundary (without includin MPI). */
  su2vector<bool>
      SolidBoundary; /*!< \brief To see if a point belong to the physical boundary (without includin MPI). */
  su2vector<bool>
      ViscousBoundary; /*!< \brief To see if a point belong to the physical boundary (without includin MPI). */
  su2vector<bool>
      PeriodicBoundary; /*!< \brief To see if a point belongs to a periodic boundary (without including MPI). */

  su2activematrix Coord; /*!< \brief vector with the coordinates of the node. */
  su2activematrix
      Coord_Old; /*!< \brief Old coordinates vector for primal solution reloading for Disc.Adj. with dynamic grid. */
  su2activematrix Coord_Sum; /*!< \brief Sum of coordinates vector for geometry smoothing. */
  su2activematrix Coord_n;   /*!< \brief Coordinates at time n for use with dynamic meshes. */
  su2activematrix Coord_n1;  /*!< \brief Coordinates at time n-1 for use with dynamic meshes. */
  su2activematrix Coord_p1;  /*!< \brief Coordinates at time n+1 for use with dynamic meshes. */

  su2activematrix GridVel;      /*!< \brief Velocity of the grid for dynamic mesh cases. */
  CVectorOfMatrix GridVel_Grad; /*!< \brief Gradient of the grid velocity for dynamic meshes. */

  su2vector<unsigned long> Parent_CV; /*!< \brief Index of the parent control volume in the agglomeration process. */
  su2vector<unsigned short> nChildren_CV; /*!< \brief Number of children in the agglomeration process. */
  vector<vector<unsigned long> >
      Children_CV; /*!< \brief Index of the children control volumes in the agglomeration process. */
  su2vector<bool> Agglomerate_Indirect; /*!< \brief This flag indicates if the indirect points can be agglomerated. */
  su2vector<bool> Agglomerate;          /*!< \brief This flag indicates if the element has been agglomerated. */

  su2vector<unsigned short> nNeighbor; /*!< \brief Number of neighbors, needed by some numerical methods. */

  /*--- Closest element on a viscous wall, and distance to it. ---*/
  su2activevector Wall_Distance;              /*!< \brief Distance to the nearest wall. */
  su2vector<int> ClosestWall_Rank;            /*!< \brief Rank of process holding the closest wall element. */
  su2vector<unsigned short> ClosestWall_Zone; /*!< \brief Zone index of closest wall element. */
  su2vector<unsigned short>
      ClosestWall_Marker; /*!< \brief Marker index of closest wall element, for given rank and zone index. */
  su2vector<unsigned long>
      ClosestWall_Elem; /*!< \brief Element index of closest wall element, for givenrank, zone and marker index. */

  su2activevector SharpEdge_Distance; /*!< \brief Distance to a sharp edge. */
  su2activevector Curvature;          /*!< \brief Value of the surface curvature (SU2_GEO). */
  su2activevector MaxLength;          /*!< \brief The maximum cell-center to cell-center length. */
  su2activevector RoughnessHeight;    /*!< \brief Roughness of the nearest wall. */

  su2matrix<int> AD_InputIndex; /*!< \brief Indices of Coord variables in the adjoint vector. */
  su2matrix<int>
      AD_OutputIndex; /*!< \brief Indices of Coord variables in the adjoint vector after having been updated. */

  /*!
   * \brief Allocate fields required by the minimal constructor.
   */
  void MinimalAllocation(unsigned long npoint);

 public:
  /*!
   * \brief "Full" constructor of the class.
   * \param[in] npoint - Number of points (dual volumes) in the problem.
   * \param[in] ndim - Number of spatial dimensions of the problem.
   * \param[in] imesh - Index of the grid allocating the points.
   * \param[in] config - Definition of the particular problem.
   */
  CPoint(unsigned long npoint, unsigned long ndim, unsigned short imesh, const CConfig* config);

  /*!
   * \brief Minimal constructor, only allocates the structures required to read and partition a mesh file.
   * \param[in] npoint - Number of points (dual volumes) in the problem.
   * \param[in] ndim - Number of spatial dimensions of the problem.
   */
  CPoint(unsigned long npoint, unsigned long ndim);

  /*!
   * \brief Default construction is not allowed.
   */
  CPoint() = delete;

  /*!
   * \brief Allocate the variables not covered by the minimal constructor.
   */
  void FullAllocation(unsigned short imesh, const CConfig* config);

  /*!
   * \brief Get the coordinates dor the control volume.
   * \param[in] iPoint - Index of the point.
   * \param[in] iDim - Number of dimensions of the problem.
   * \return Coordinate that correspond with <i>iDim</i>.
   */
  inline su2double GetCoord(unsigned long iPoint, unsigned long iDim) const { return Coord(iPoint, iDim); }

  /*!
   * \brief Get the coordinates of the control volume.
   * \param[in] iPoint - Index of the point.
   * \return pointer to the coordinate of the point.
   */
  inline su2double* GetCoord(unsigned long iPoint) { return Coord[iPoint]; }

  /*!
   * \brief Get the entire matrix of coordinates of the control volumes.
   */
  inline const su2activematrix& GetCoord() const { return Coord; }

  /*!
   * \brief Set the coordinates for the control volume.
   * \param[in] iPoint - Index of the point.
   * \param[in] iDim - Position to store the coordinate.
   * \param[in] coord - Coordinate for iDim.
   */
  inline void SetCoord(unsigned long iPoint, unsigned long iDim, su2double coord) { Coord(iPoint, iDim) = coord; }

  /*!
   * \brief Set the coordinates for the control volume.
   * \param[in] iPoint - Index of the point.
   * \param[in] iDim - Position to store the coordinate.
   * \param[in] coord - Coordinate for iDim.
   */
  inline void AddCoord(unsigned long iPoint, unsigned long iDim, su2double coord) { Coord(iPoint, iDim) += coord; }

  /*!
   * \brief Set the point coordinates.
   * \param[in] iPoint - Index of the point.
   * \param[in] coord - Coordinate of the point.
   */
  inline void SetCoord(unsigned long iPoint, const su2double* coord) {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Coord(iPoint, iDim) = coord[iDim];
  }

  /*!
   * \brief Set the elements that are connected to each point.
   * \param[in] elemsMatrix - List of lists with the neighbor points connected to each point.
   */
  void SetElems(const vector<vector<long> >& elemsMatrix);

  /*!
   * \brief Reset the elements of a control volume.
   */
  inline void ResetElems() { Elem = CCompressedSparsePatternL(); }

  /*!
   * \brief Get the number of elements that compose the control volume.
   * \param[in] iPoint - Index of the point.
   * \return Number of elements that compose the control volume.
   */
  inline unsigned short GetnElem(unsigned long iPoint) const { return Elem.getNumNonZeros(iPoint); }

  /*!
   * \brief Get all the elements that compose the control volume.
   * \param[in] iPoint - Index of the point.
   * \param[in] nelem - Position where the element is stored.
   * \return Index of the element.
   */
  inline unsigned long GetElem(unsigned long iPoint, unsigned long nelem) const {
    return Elem.getInnerIdx(iPoint, nelem);
  }

  /*!
   * \brief Get inner iterator to loop over neighbor elements.
   */
  inline CCompressedSparsePatternL::CInnerIter GetElems(unsigned long iPoint) const {
    return Elem.getInnerIter(iPoint);
  }

  /*!
   * \brief Set the points that compose the control volume.
   * \param[in] pointsMatrix - List of lists with the neighbor points connected to each point.
   */
  void SetPoints(const vector<vector<unsigned long> >& pointsMatrix);

  /*!
   * \brief Get the entire point adjacency information in compressed format (CSR).
   */
  const CCompressedSparsePatternUL& GetPoints() const { return Point; }

  /*!
   * \brief Reset the points that compose the control volume.
   */
  inline void ResetPoints() {
    Point = CCompressedSparsePatternUL();
    Edge = CCompressedSparsePatternL();
  }

  /*!
   * \brief Get the number of points that compose the control volume.
   * \param[in] iPoint - Index of the point.
   * \return Number of points that compose the control volume.
   */
  inline unsigned short GetnPoint(unsigned long iPoint) const { return Point.getNumNonZeros(iPoint); }

  /*!
   * \brief Get all the points that compose the control volume.
   * \param[in] iPoint - Index of the point.
   * \param[in] point - Position where the point is stored.
   * \return Index of the point.
   */
  inline unsigned long GetPoint(unsigned long iPoint, unsigned long npoint) const {
    return Point.getInnerIdx(iPoint, npoint);
  }

  /*!
   * \brief Get inner iterator to loop over neighbor points.
   */
  inline CCompressedSparsePatternUL::CInnerIter GetPoints(unsigned long iPoint) const {
    return Point.getInnerIter(iPoint);
  }

  /*!
   * \brief Set the edges that compose the control volume.
   * \param[in] iPoint - Index of the point.
   * \param[in] iedge - Edge to be added.
   * \param[in] nedge - Position in which is going to be stored the edge for each control volume.
   */
  inline void SetEdge(unsigned long iPoint, long iedge, unsigned long nedge) {
    Edge.getInnerIdx(iPoint, nedge) = iedge;
  }

  /*!
   * \brief Get all the edges that compose the control volume.
   * \param[in] iPoint - Index of the point.
   * \param[in] nedge - Position where the edge is stored.
   * \return Index of the edge.
   */
  inline long GetEdge(unsigned long iPoint, unsigned long nedge) const { return Edge.getInnerIdx(iPoint, nedge); }

  /*!
   * \brief Get inner iterator to loop over neighbor edges.
   */
  inline CCompressedSparsePatternL::CInnerIter GetEdges(unsigned long iPoint) const {
    return Edge.getInnerIter(iPoint);
  }

  /*!
   * \brief Set the boundary vertex that compose the control volume.
   * \param[in] iPoint - Index of the point.
   * \param[in] iVertex - Vertex to be added.
   * \param[in] iMarker - Marker of the vertex to be added (position where is going to be stored).
   */
  inline void SetVertex(unsigned long iPoint, long iVertex, unsigned long iMarker) {
    if (Boundary(iPoint)) Vertex[iPoint][iMarker] = iVertex;
  }

  /*!
   * \brief Get the vertex that compose the control volume for a marker.
   * \param[in] iPoint - Index of the point.
   * \param[in] iMarker - Position where the vertex is stored.
   * \return Index of the vertex.
   */
  inline long GetVertex(unsigned long iPoint, unsigned long iMarker) const {
    if (Boundary(iPoint))
      return Vertex[iPoint][iMarker];
    else
      return -1;
  }

  /*!
   * \brief Set if a point belong to the boundary.
   * \param[in] iPoint - Index of the point.
   * \note It also create the structure to store the vertex.
   * \param[in] nMarker - Max number of marker.
   */
  inline void SetBoundary(unsigned long iPoint, unsigned short nMarker) {
    if (!Boundary(iPoint)) Vertex[iPoint].resize(nMarker, -1);
    Boundary(iPoint) = true;
  }

  /*!
   * \brief Reset the boundary of a control volume.
   * \param[in] iPoint - Index of the point.
   */
  inline void ResetBoundary(unsigned long iPoint) {
    Vertex[iPoint].clear();
    Boundary(iPoint) = false;
  }

  /*!
   * \brief Mark the point as boundary.
   * \param[in] iPoint - Index of the point.
   * \param[in] boundary - <code>TRUE</code> if the point belong to the boundary; otherwise <code>FALSE</code>.
   */
  inline void SetBoundary(unsigned long iPoint, bool boundary) { Boundary(iPoint) = boundary; }

  /*!
   * \brief Provides information about if a point belong to the boundaries.
   * \param[in] iPoint - Index of the point.
   * \return <code>TRUE</code> if the point belong to the boundary; otherwise <code>FALSE</code>.
   */
  inline bool GetBoundary(unsigned long iPoint) const { return Boundary(iPoint); }

  /*!
   * \brief Set if a point belong to the boundary.
   * \param[in] iPoint - Index of the point.
   * \param[in] boundary - <code>TRUE</code> if the point belong to the physical boundary; otherwise <code>FALSE</code>.
   */
  inline void SetPhysicalBoundary(unsigned long iPoint, bool boundary) { PhysicalBoundary(iPoint) = boundary; }

  /*!
   * \brief Provides information about if a point belong to the physical boundaries (without MPI).
   * \param[in] iPoint - Index of the point.
   * \return <code>TRUE</code> if the point belong to the boundary; otherwise <code>FALSE</code>.
   */
  inline bool GetPhysicalBoundary(unsigned long iPoint) const { return PhysicalBoundary(iPoint); }

  /*!
   * \brief Set if a point belong to the boundary.
   * \param[in] iPoint - Index of the point.
   * \param[in] boundary - <code>TRUE</code> if the point belong to the physical boundary; otherwise <code>FALSE</code>.
   */
  inline void SetSolidBoundary(unsigned long iPoint, bool boundary) { SolidBoundary(iPoint) = boundary; }

  /*!
   * \brief Provides information about if a point belong to the physical boundaries (without MPI).
   * \param[in] iPoint - Index of the point.
   * \return <code>TRUE</code> if the point belong to the boundary; otherwise <code>FALSE</code>.
   */
  inline bool GetSolidBoundary(unsigned long iPoint) const { return SolidBoundary(iPoint); }

  /*!
   * \brief Set if a point belong to the boundary.
   * \param[in] iPoint - Index of the point.
   * \param[in] boundary - <code>TRUE</code> if the point belong to the physical boundary; otherwise <code>FALSE</code>.
   */
  inline void SetViscousBoundary(unsigned long iPoint, bool boundary) { ViscousBoundary(iPoint) = boundary; }

  /*!
   * \brief Provides information about if a point belong to the physical boundaries (without MPI).
   * \param[in] iPoint - Index of the point.
   * \return <code>TRUE</code> if the point belong to the boundary; otherwise <code>FALSE</code>.
   */
  inline bool GetViscousBoundary(unsigned long iPoint) const { return ViscousBoundary(iPoint); }

  /*!
   * \brief Set if a point belongs to a periodic boundary.
   * \param[in] iPoint - Index of the point.
   * \param[in] boundary - <code>TRUE</code> if the point belongs to a periodic boundary; otherwise <code>FALSE</code>.
   */
  inline void SetPeriodicBoundary(unsigned long iPoint, bool boundary) { PeriodicBoundary(iPoint) = boundary; }

  /*!
   * \brief Provides information about if a point belongs to a periodic boundary (without MPI).
   * \param[in] iPoint - Index of the point.
   * \return <code>TRUE</code> if the point belongs to a periodic boundary; otherwise <code>FALSE</code>.
   */
  inline bool GetPeriodicBoundary(unsigned long iPoint) const { return PeriodicBoundary(iPoint); }

  /*!
   * \brief For parallel computation, its indicates if a point must be computed or not.
   * \param[in] iPoint - Index of the point.
   * \param[in] domain - <code>TRUE</code> if the point belong to the domain; otherwise <code>FALSE</code>.
   */
  inline void SetDomain(unsigned long iPoint, bool domain) { Domain(iPoint) = domain; }

  /*!
   * \brief For parallel computation, its indicates if a point must be computed or not.
   * \param[in] iPoint - Index of the point.
   * \return <code>TRUE</code> if the node belong to the physical domain; otherwise <code>FALSE</code>.
   */
  inline bool GetDomain(unsigned long iPoint) const { return Domain(iPoint); }

  /*!
   * \brief Set a color to the point that comes from the grid partitioning.
   * \param[in] iPoint - Index of the point.
   * \note Each domain has a different color.
   * \param[in] color - Color of the point.
   */
  inline void SetColor(unsigned long iPoint, unsigned long color) { Color(iPoint) = color; }

  /*!
   * \brief Get the color of a point, the color indicates to which subdomain the point belong to.
   * \param[in] iPoint - Index of the point.
   * \return Color of the point.
   */
  inline unsigned long GetColor(unsigned long iPoint) const { return Color(iPoint); }

  /*!
   * \brief Set the global index in a parallel computation.
   * \param[in] iPoint - Index of the point.
   * \return Global index in a parallel computation.
   */
  inline void SetGlobalIndex(unsigned long iPoint, unsigned long globalindex) { GlobalIndex(iPoint) = globalindex; }

  /*!
   * \brief Get the global index in a parallel computation.
   * \param[in] iPoint - Index of the point.
   * \return Global index in a parallel computation.
   */
  inline unsigned long GetGlobalIndex(unsigned long iPoint) const { return GlobalIndex(iPoint); }

  /*!
   * \brief Set the number of neighbor (artificial dissipation).
   * \param[in] iPoint - Index of the point.
   * \param[in] nneighbor - Number of neighbors.
   */
  inline void SetnNeighbor(unsigned long iPoint, unsigned short nneighbor) { nNeighbor(iPoint) = nneighbor; }

  /*!
   * \brief Get the number of neighbor of a point.
   * \param[in] iPoint - Index of the point.
   * \return Number of neighbors.
   */
  inline unsigned short GetnNeighbor(unsigned long iPoint) const { return nNeighbor(iPoint); }

  /*!
   * \brief Set the value of the distance to the nearest wall.
   * \param[in] iPoint - Index of the point.
   * \param[in] distance - Value of the distance.
   * \param[in] rankID - Rank of process holding the closest wall element.
   * \param[in] zoneID - Zone index of closest wall element.
   * \param[in] markerID - Marker index of closest wall element.
   * \param[in] elemID - Element index of closest wall element.
   */
  inline void SetWall_Distance(unsigned long iPoint, su2double distance, int rankID, unsigned short zoneID,
                               unsigned short markerID, unsigned long elemID) {
    Wall_Distance(iPoint) = distance;
    ClosestWall_Rank(iPoint) = rankID;
    ClosestWall_Zone(iPoint) = zoneID;
    ClosestWall_Marker(iPoint) = markerID;
    ClosestWall_Elem(iPoint) = elemID;
  }
  inline void SetWall_Distance(unsigned long iPoint, su2double distance) { Wall_Distance(iPoint) = distance; }

  /*!
   * \brief Get the value of the distance to the nearest wall.
   * \param[in] iPoint - Index of the point.
   * \return Value of the distance to the nearest wall.
   */
  inline su2double& GetWall_Distance(unsigned long iPoint) { return Wall_Distance(iPoint); }
  inline const su2double& GetWall_Distance(unsigned long iPoint) const { return Wall_Distance(iPoint); }

  /*!
   * \brief Set the value of the distance to the nearest wall.
   * \param[in] iPoint - Index of the point.
   * \param[in] distance - Value of the distance.
   */
  inline void SetRoughnessHeight(unsigned long iPoint, su2double roughheight) { RoughnessHeight(iPoint) = roughheight; }

  /*!
   * \brief Get the value of the distance to the nearest wall.
   * \param[in] iPoint - Index of the point.
   * \return Value of the distance to the nearest wall.
   */
  inline su2double GetRoughnessHeight(unsigned long iPoint) const { return RoughnessHeight(iPoint); }

  /*!
   * \brief Set the value of the distance to a sharp edge.
   * \param[in] iPoint - Index of the point.
   * \param[in] distance - Value of the distance.
   */
  inline void SetSharpEdge_Distance(unsigned long iPoint, su2double distance) { SharpEdge_Distance(iPoint) = distance; }

  /*!
   * \brief Get the value of the distance to a sharp edge
   * \param[in] iPoint - Index of the point.
   * \return Value of the distance to the nearest wall.
   */
  inline su2double& GetSharpEdge_Distance(unsigned long iPoint) { return SharpEdge_Distance(iPoint); }
  inline const su2double& GetSharpEdge_Distance(unsigned long iPoint) const { return SharpEdge_Distance(iPoint); }

  /*!
   * \brief Set the value of the curvature at a surface node.
   * \param[in] iPoint - Index of the point.
   * \param[in] curvature - Value of the curvature.
   */
  inline void SetCurvature(unsigned long iPoint, su2double curvature) { Curvature(iPoint) = curvature; }

  /*!
   * \brief Get the value of the curvature at a surface node.
   * \param[in] iPoint - Index of the point.
   * \return Value of the curvature.
   */
  inline su2double GetCurvature(unsigned long iPoint) const { return Curvature(iPoint); }

  /*!
   * \brief Set the max cell-center to cell-center length.
   * \param[in] iPoint - Index of the point.
   * \param[in] max_length - Value of the max length
   */
  inline void SetMaxLength(unsigned long iPoint, su2double max_length) { MaxLength(iPoint) = max_length; }

  /*!
   * \brief Get the maximum cell-center to cell-center length.
   * \param[in] iPoint - Index of the point.
   * \return The maximum cell-center to cell-center length.
   */
  inline su2double GetMaxLength(unsigned long iPoint) const { return MaxLength(iPoint); }

  /*!
   * \brief Get area or volume of the control volume.
   * \param[in] iPoint - Index of the point.
   * \return Area or volume of the control volume.
   */
  inline su2double& GetVolume(unsigned long iPoint) { return Volume(iPoint); }
  inline const su2double& GetVolume(unsigned long iPoint) const { return Volume(iPoint); }

  /*!
   * \brief Set the volume of the control volume.
   * \param[in] iPoint - Index of the point.
   * \param[in] volume - Value of the volume.
   */
  inline void SetVolume(unsigned long iPoint, su2double volume) { Volume(iPoint) = volume; }

  /*!
   * \brief Adds some area or volume of the CV.
   * \param[in] iPoint - Index of the point.
   * \param[in] volume - Local volume to be added to the total one.
   */
  inline void AddVolume(unsigned long iPoint, su2double volume) { Volume(iPoint) += volume; }

  /*!
   * \brief Get the missing component of area or volume for a control volume on a periodic marker.
   * \param[in] iPoint - Index of the point.
   * \return Periodic component of area or volume for a control volume on a periodic marker.
   */
  inline su2double& GetPeriodicVolume(unsigned long iPoint) { return Periodic_Volume(iPoint); }
  inline const su2double& GetPeriodicVolume(unsigned long iPoint) const { return Periodic_Volume(iPoint); }

  /*!
   * \brief Set the missing component of area or volume for a control volume on a periodic marker.
   * \param[in] iPoint - Index of the point.
   * \param[in] volume - Value of the volume from the missing components of the CV on the periodic marker.
   */
  inline void SetPeriodicVolume(unsigned long iPoint, su2double volume) { Periodic_Volume(iPoint) = volume; }

  /*!
   * \brief Get the volume of the control volume at time n.
   * \param[in] iPoint - Index of the point.
   * \return Volume of the control volume at time n
   */
  inline su2double GetVolume_n(unsigned long iPoint) const { return Volume_n(iPoint); }

  /*!
   * \brief Get the volume of the control volume at time n+1.
   * \param[in] iPoint - Index of the point.
   * \return Volume of the control volume at time n+1
   */
  inline su2double GetVolume_nM1(unsigned long iPoint) const { return Volume_nM1(iPoint); }

  /*!
   * \brief Set the volume of the control volume at time n.
   */
  void SetVolume_n();

  /*!
   * \brief Set the volume of the control volume at time n-1.
   */
  void SetVolume_nM1();

  /*!
   * \brief Set the volume of the control volume at time n using n-1.
   */
  void SetVolume_n_from_OldnM1();

  /*!
   * \brief Set the volume of the control volume at current time using time n.
   */
  void SetVolume_from_Oldn();

  /*!
   * \brief Set the Volume to Volume_Old.
   */
  void SetVolume_Old();

  /*!
   * \brief Set the Volume_n to Volume_n_Old.
   */
  void SetVolume_n_Old();

  /*!
   * \brief Set the Volume_nM1 to Volume_nM1_Old.
   */
  void SetVolume_nM1_Old();

  /*!
   * \brief Set the parent control volume of an agglomerated control volume.
   * \param[in] iPoint - Index of the point.
   * \param[in] parent_CV - Index of the parent control volume.
   */
  inline void SetParent_CV(unsigned long iPoint, unsigned long parent_CV) {
    Parent_CV(iPoint) = parent_CV;
    Agglomerate(iPoint) = true;
  }

  /*!
   * \brief Set the children control volumes of an agglomerated control volume.
   * \param[in] iPoint - Index of the point.
   * \param[in] nchildren_CV - Number of children.
   * \param[in] children_CV - Index of the children control volume.
   */
  inline void SetChildren_CV(unsigned long iPoint, unsigned long nchildren_CV, unsigned long children_CV) {
    Children_CV[iPoint].resize(nchildren_CV + 1);
    Children_CV[iPoint][nchildren_CV] = children_CV;
  }

  /*!
   * \brief Get the parent control volume of an agglomerated control volume.
   * \param[in] iPoint - Index of the point.
   * \return Index of the parent control volume.
   */
  inline unsigned long GetParent_CV(unsigned long iPoint) const { return Parent_CV(iPoint); }

  /*!
   * \brief Get the children control volume of an agglomerated control volume.
   * \param[in] iPoint - Index of the point.
   * \param[in] nchildren_CV - Number of the children.
   * \return Index of the parent control volume.
   */
  inline unsigned long GetChildren_CV(unsigned long iPoint, unsigned short nchildren_CV) const {
    return Children_CV[iPoint][nchildren_CV];
  }

  /*!
   * \brief Get information about if a control volume has been agglomerated.
   * \param[in] iPoint - Index of the point.
   * \return <code>TRUE</code> if the point has been agglomerated; otherwise <code>FALSE</code>.
   */
  inline bool GetAgglomerate(unsigned long iPoint) const { return Agglomerate(iPoint); }

  /*!
   * \brief Get information about if the indirect neighbors can be agglomerated.
   * \param[in] iPoint - Index of the point.
   * \return <code>TRUE</code> if the indirect neigbors can be agglomerated; otherwise <code>FALSE</code>.
   */
  inline bool GetAgglomerate_Indirect(unsigned long iPoint) const { return Agglomerate_Indirect(iPoint); }

  /*!
   * \brief Set information about if the indirect neighbors can be agglomerated.
   * \param[in] iPoint - Index of the point.
   * \param[in] agglomerate - The indirect neigbors can be agglomerated.
   */
  inline void SetAgglomerate_Indirect(unsigned long iPoint, bool agglomerate) {
    Agglomerate_Indirect(iPoint) = agglomerate;
  };

  /*!
   * \brief Get the number of children of an agglomerated control volume.
   * \param[in] iPoint - Index of the point.
   * \return Number of children control volume.
   */
  inline unsigned short GetnChildren_CV(unsigned long iPoint) const { return nChildren_CV(iPoint); }

  /*!
   * \brief Set the number of children of an agglomerated control volume.
   * \param[in] iPoint - Index of the point.
   * \param[in] nchildren_CV - Number of children of the control volume.
   */
  inline void SetnChildren_CV(unsigned long iPoint, unsigned short nchildren_CV) {
    nChildren_CV(iPoint) = nchildren_CV;
  }

  /*!
   * \brief Get the coordinates of the control volume at time n.
   * \param[in] iPoint - Index of the point.
   * \return Coordinates of the control volume at time n.
   */
  inline su2double* GetCoord_n(unsigned long iPoint) { return Coord_n[iPoint]; }

  /*!
   * \brief Get the coordinates of the control volume at time n-1.
   * \param[in] iPoint - Index of the point.
   * \return Volume of the control volume at time n-1
   */
  inline su2double* GetCoord_n1(unsigned long iPoint) { return Coord_n1[iPoint]; }

  /*!
   * \brief Get the coordinates of the control volume at time n+1.
   * \param[in] iPoint - Index of the point.
   * \return Volume of the control volume at time n+1
   */
  inline su2double* GetCoord_p1(unsigned long iPoint) { return Coord_p1[iPoint]; }

  /*!
   * \brief Set the coordinates of the control volume at time n to the ones in <i>Coord</i>.
   */
  void SetCoord_n();

  /*!
   * \brief Set the coordinates of the control volume at time n-1 to the ones in <i>Coord_n</i>.
   */
  void SetCoord_n1();

  /*!
   * \brief Set the coordinates of the control volume at time n, for restart cases.
   * \param[in] iPoint - Index of the point.
   * \param[in] coord - Value of the grid coordinates at time n.
   */
  inline void SetCoord_n(unsigned long iPoint, const su2double* coord) {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Coord_n(iPoint, iDim) = coord[iDim];
  }

  /*!
   * \brief Set the coordinates of the control volume at time n-1, for restart cases.
   * \param[in] iPoint - Index of the point.
   * \param[in] coord - Value of the grid coordinates at time n-1.
   */
  inline void SetCoord_n1(unsigned long iPoint, const su2double* coord) {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Coord_n1(iPoint, iDim) = coord[iDim];
  }

  /*!
   * \brief Set the coordinates of the control volume at time n+1.
   * \param[in] iPoint - Index of the point.
   * \param[in] coord - Value of the grid coordinates at time n+1.
   */
  inline void SetCoord_p1(unsigned long iPoint, const su2double* coord) {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Coord_p1(iPoint, iDim) = coord[iDim];
  }

  /*!
   * \brief Get the value of the old coordinates for implicit smoothing.
   * \param[in] iPoint - Index of the point.
   * \return Old coordinates at a point.
   */
  inline su2double* GetCoord_Old(unsigned long iPoint) { return Coord_Old[iPoint]; }

  /*!
   * \brief Set the value of the vector <i>Coord_Old</i> for implicit smoothing.
   * \param[in] iPoint - Index of the point.
   * \param[in] coord_old - Value of the coordinates.
   */
  inline void SetCoord_Old(unsigned long iPoint, const su2double* coord_old) {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Coord_Old(iPoint, iDim) = coord_old[iDim];
  }

  /*!
   * \brief Set the value of the vector <i>Coord_Old</i> to <i>Coord</i>.
   */
  void SetCoord_Old();

  /*!
   * \brief Get the value of the summed coordinates for implicit smoothing.
   * \param[in] iPoint - Index of the point.
   * \return Sum of coordinates at a point.
   */
  inline su2double* GetCoord_Sum(unsigned long iPoint) { return Coord_Sum[iPoint]; }

  /*!
   * \brief Add the value of the coordinates to the <i>Coord_Sum</i> vector for implicit smoothing.
   * \param[in] iPoint - Index of the point.
   * \param[in] coord_sum - Value of the coordinates to add.
   */
  inline void AddCoord_Sum(unsigned long iPoint, const su2double* coord_sum) {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Coord_Sum(iPoint, iDim) += coord_sum[iDim];
  }

  /*!
   * \brief Initialize the vector <i>Coord_Sum</i>.
   */
  void SetCoord_SumZero();

  /*!
   * \brief Get the value of the grid velocity at the point.
   * \param[in] iPoint - Index of the point.
   * \return Grid velocity at the point.
   */
  inline su2double* GetGridVel(unsigned long iPoint) { return GridVel[iPoint]; }

  /*!
   * \brief Get the grid velocity matrix for the entire domain.
   */
  inline const su2activematrix& GetGridVel() const { return GridVel; }

  /*!
   * \brief Set the value of the grid velocity at the point.
   * \param[in] iPoint - Index of the point.
   * \param[in] iDim - Index of the coordinate.
   * \param[in] gridvel - Value of the grid velocity.
   */
  inline void SetGridVel(unsigned long iPoint, unsigned long iDim, su2double gridvel) {
    GridVel(iPoint, iDim) = gridvel;
  }

  /*!
   * \brief Set the value of the grid velocity at the point.
   * \param[in] iPoint - Index of the point.
   * \param[in] gridvel - Value of the grid velocity.
   */
  inline void SetGridVel(unsigned long iPoint, const su2double* gridvel) {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) GridVel(iPoint, iDim) = gridvel[iDim];
  }

  /*!
   * \brief Get the grid velocity gradients for the entire domain.
   */
  inline CVectorOfMatrix& GetGridVel_Grad() { return GridVel_Grad; }

  /*!
   * \brief Get the value of the grid velocity gradient at the point.
   * \param[in] iPoint - Index of the point.
   * \return Grid velocity gradient at the point.
   */
  inline CMatrixView<const su2double> GetGridVel_Grad(unsigned long iPoint) const { return GridVel_Grad[iPoint]; }

  /*!
   * \brief Set the adjoint values of the (geometric) coordinates.
   * \param[in] iPoint - Index of the point.
   * \param[in] adj_sol - Adjoint values of the Coord variables.
   */
  inline void SetAdjointSolution(unsigned long iPoint, const su2double* adj_sol) {
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      AD::SetDerivative(AD_OutputIndex(iPoint, iDim), SU2_TYPE::GetValue(adj_sol[iDim]));
  }

  /*!
   * \brief Get the adjoint values of the (geometric) coordinates.
   * \param[in] iPoint - Index of the point.
   * \param[in] iDim - Dimension.
   */
  inline su2double GetAdjointSolution(unsigned long iPoint, unsigned long iDim) const {
    return AD::GetDerivative(AD_InputIndex(iPoint, iDim));
  }

  /*!
   * \brief Register coordinates of a point.
   * \param[in] iPoint - Index of the point.
   * \param[in] input - Register as input or output.
   */
  inline void RegisterCoordinates(unsigned long iPoint, bool input) {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) {
      if (input) {
        AD::RegisterInput(Coord(iPoint, iDim));
        AD::SetIndex(AD_InputIndex(iPoint, iDim), Coord(iPoint, iDim));
      } else {
        AD::RegisterOutput(Coord(iPoint, iDim));
        AD::SetIndex(AD_OutputIndex(iPoint, iDim), Coord(iPoint, iDim));
      }
    }
  }

  /*!
   * \brief Set wall roughnesses according to stored closest wall information.
   * \param[in] roughness - Mapping [rank][zone][marker] -> roughness
   */
  template <typename Roughness_type>
  void SetWallRoughness(Roughness_type const& roughness) {
    for (unsigned long iPoint = 0; iPoint < GlobalIndex.size(); ++iPoint) {
      auto rankID = ClosestWall_Rank[iPoint];
      auto zoneID = ClosestWall_Zone[iPoint];
      auto markerID = ClosestWall_Marker[iPoint];
      if (rankID >= 0) {
        SetRoughnessHeight(iPoint, roughness[rankID][zoneID][markerID]);
      }
    }
  }
};
