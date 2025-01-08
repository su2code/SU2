/*!
 * \file CGeometry.hpp
 * \brief Headers of the main subroutines for creating the geometrical structure.
 *        The subroutines and functions are in the <i>CGeometry.cpp</i> file.
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

#include <limits>
#include "../parallelization/mpi_structure.hpp"

#ifdef HAVE_METIS
#include "metis.h"
#endif
#ifdef HAVE_PARMETIS
extern "C" {
#include "parmetis.h"
}
#endif

#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <memory>
#include <unordered_map>

#include "primal_grid/CPrimalGrid.hpp"
#include "dual_grid/CDualGrid.hpp"
#include "dual_grid/CPoint.hpp"
#include "dual_grid/CEdge.hpp"
#include "dual_grid/CVertex.hpp"
#include "dual_grid/CTurboVertex.hpp"

#include "../CConfig.hpp"
#include "../fem/geometry_structure_fem_part.hpp"
#include "../toolboxes/graph_toolbox.hpp"
#include "../adt/CADTElemClass.hpp"

using namespace std;

/*!
 * \class CGeometry
 * \brief Parent class for defining the geometry of the problem (complete geometry,
 *        multigrid agglomerated geometry, only boundary geometry, etc..)
 * \author F. Palacios
 */
class CGeometry {
 protected:
  enum : size_t { OMP_MIN_SIZE = 32 }; /*!< \brief Chunk size for small loops. */
  enum : size_t { MAXNDIM = 3 };

  const int size{SINGLE_NODE}; /*!< \brief MPI Size. */
  const int rank{MASTER_NODE}; /*!< \brief MPI Rank. */

  unsigned long nPoint{0}, /*!< \brief Number of points of the mesh. */
      nPointDomain{0},     /*!< \brief Number of real points of the mesh. */
      nPointGhost{0},      /*!< \brief Number of ghost points of the mesh. */
      Global_nPoint{0},    /*!< \brief Total number of nodes in a simulation across all processors (including halos). */
      Global_nPointDomain{
          0},          /*!< \brief Total number of nodes in a simulation across all processors (excluding halos). */
      nElem{0},        /*!< \brief Number of elements of the mesh. */
      Global_nElem{0}, /*!< \brief Total number of elements in a simulation across all processors (all types). */
      Global_nElemDomain{
          0},        /*!< \brief Total number of elements in a simulation across all processors (excluding halos). */
      nEdge{0},      /*!< \brief Number of edges of the mesh. */
      nFace{0},      /*!< \brief Number of faces of the mesh. */
      nelem_edge{0}, /*!< \brief Number of edges in the mesh. */
      Global_nelem_edge{0},       /*!< \brief Total number of edges in the mesh across all processors. */
      nelem_triangle{0},          /*!< \brief Number of triangles in the mesh. */
      Global_nelem_triangle{0},   /*!< \brief Total number of triangles in the mesh across all processors. */
      nelem_quad{0},              /*!< \brief Number of quadrangles in the mesh. */
      Global_nelem_quad{0},       /*!< \brief Total number of quadrangles in the mesh across all processors. */
      nelem_tetra{0},             /*!< \brief Number of tetrahedra in the mesh. */
      Global_nelem_tetra{0},      /*!< \brief Total number of tetrahedra in the mesh across all processors. */
      nelem_hexa{0},              /*!< \brief Number of hexahedra in the mesh. */
      Global_nelem_hexa{0},       /*!< \brief Total number of hexahedra in the mesh across all processors. */
      nelem_prism{0},             /*!< \brief Number of prisms in the mesh. */
      Global_nelem_prism{0},      /*!< \brief Total number of prisms in the mesh across all processors. */
      nelem_pyramid{0},           /*!< \brief Number of pyramids in the mesh. */
      Global_nelem_pyramid{0},    /*!< \brief Total number of pyramids in the mesh across all processors. */
      nelem_edge_bound{0},        /*!< \brief Number of edges on the mesh boundaries. */
      Global_nelem_edge_bound{0}, /*!< \brief Total number of edges on the mesh boundaries across all processors. */
      nelem_triangle_bound{0},    /*!< \brief Number of triangles on the mesh boundaries. */
      Global_nelem_triangle_bound{
          0},                     /*!< \brief Total number of triangles on the mesh boundaries across all processors. */
      nelem_quad_bound{0},        /*!< \brief Number of quads on the mesh boundaries. */
      Global_nelem_quad_bound{0}, /*!< \brief Total number of quads on the mesh boundaries across all processors. */
      nNonconvexElements{0};      /*!< \brief Number of nonconvex elements in the mesh. */

  unsigned short nDim{0};       /*!< \brief Number of dimension of the problem. */
  unsigned short nZone{0};      /*!< \brief Number of zones in the problem. */
  unsigned short nMarker{0};    /*!< \brief Number of different markers of the mesh. */
  unsigned short nCommLevel{0}; /*!< \brief Number of non-blocking communication levels. */

  unsigned short MGLevel{0};        /*!< \brief The mesh level index for the current geometry container. */
  unsigned long Max_GlobalPoint{0}; /*!< \brief Greater global point in the domain local structure. */

  /*--- Boundary information. ---*/

  short* Marker_All_SendRecv{nullptr}; /*!< \brief MPI Marker. */
  su2double** CustomBoundaryTemperature{nullptr};
  su2double** CustomBoundaryHeatFlux{nullptr};

  /*--- Create vectors and distribute the values among the different planes queues ---*/

  vector<vector<su2double>>
      Xcoord_plane; /*!< \brief Vector containing x coordinates of new points appearing on a single plane */
  vector<vector<su2double>>
      Ycoord_plane; /*!< \brief Vector containing y coordinates of new points appearing on a single plane */
  vector<vector<su2double>>
      Zcoord_plane; /*!< \brief Vector containing z coordinates of new points appearing on a single plane */
  vector<vector<su2double>> FaceArea_plane;   /*!< \brief Vector containing area/volume associated with  new points
                                                 appearing on a single plane */
  vector<vector<unsigned long>> Plane_points; /*!< \brief Vector containing points appearing on a single plane */

  vector<su2double> XCoordList; /*!< \brief Vector containing points appearing on a single plane */

#if defined(HAVE_MPI) && defined(HAVE_PARMETIS)
  vector<vector<unsigned long>>
      adj_nodes; /*!< \brief Vector of vectors holding each node's adjacency during preparation for ParMETIS. */
  vector<idx_t> adjacency; /*!< \brief Local adjacency array to be input into ParMETIS for partitioning (idx_t is a
                              ParMETIS type defined in their headers). */
  vector<idx_t> xadj; /*!< \brief Index array that points to the start of each node's adjacency in CSR format (needed to
                         interpret the adjacency array).  */
#endif

  /*--- Turbomachinery variables ---*/

  unsigned short* nSpanWiseSections{
      nullptr}; /*!< \brief Number of Span wise section for each turbo marker, indexed by inflow/outflow */
  unsigned short* nSpanSectionsByMarker{nullptr}; /*!< \brief Number of Span wise section for each turbo marker, indexed
                                                     by marker.  Needed for deallocation.*/
  unsigned short nTurboPerf{0};                   /*!< \brief Number of Span wise section for each turbo marker. */
  su2double** SpanWiseValue{nullptr};             /*!< \brief Span wise values for each turbo marker. */
  long** nVertexSpan{nullptr};             /*!< \brief number of vertexes for span wise section for each marker.  */
  unsigned long** nTotVertexSpan{nullptr}; /*!< \brief number of vertexes at each span wise section for each marker.  */
  unsigned long nVertexSpanMax[3] = {
      0}; /*!< \brief max number of vertexes for each span section for each marker flag.  */
  su2double*** AverageTurboNormal{nullptr}; /*!< \brief Average boundary normal at each span wise section for each
                                               marker in the turbomachinery frame of reference.*/
  su2double*** AverageNormal{nullptr}; /*!< \brief Average boundary normal at each span wise section for each marker.*/
  su2double*** AverageGridVel{
      nullptr}; /*!< \brief Average boundary grid velocity at each span wise section for each marker.*/
  su2double** AverageTangGridVel{
      nullptr}; /*!< \brief Average tangential rotational speed at each span wise section for each marker.*/
  su2double** SpanArea{nullptr};        /*!< \brief Area at each span wise section for each marker.*/
  su2double** MaxAngularCoord{nullptr}; /*!< \brief Max angular pitch at each span wise section for each marker.*/
  su2double** MinAngularCoord{nullptr}; /*!< \brief Max angular pitch at each span wise section for each marker.*/
  su2double** MinRelAngularCoord{
      nullptr};                     /*!< \brief Min relative angular coord at each span wise section for each marker.*/
  su2double** TurboRadius{nullptr}; /*!< \brief Radius at each span wise section for each marker.*/
  su2double** TangGridVelIn{nullptr};
  su2double** TangGridVelOut{nullptr}; /*!< \brief Average tangential rotational speed at each span wise section for
                                          each turbomachinery marker.*/
  su2double** SpanAreaIn{nullptr};
  su2double** SpanAreaOut{nullptr}; /*!< \brief Area at each span wise section for each turbomachinery marker.*/
  su2double** TurboRadiusIn{nullptr};
  su2double** TurboRadiusOut{nullptr}; /*!< \brief Radius at each span wise section for each turbomachinery marker*/

  /*--- Sparsity patterns associated with the geometry. ---*/

  CCompressedSparsePatternUL finiteVolumeCSRFill0, /*!< \brief 0-fill FVM sparsity. */
      finiteVolumeCSRFillN,                        /*!< \brief N-fill FVM sparsity (e.g. for ILUn preconditioner). */
      finiteElementCSRFill0,                       /*!< \brief 0-fill FEM sparsity. */
      finiteElementCSRFillN;                       /*!< \brief N-fill FEM sparsity (e.g. for ILUn preconditioner). */

  CEdgeToNonZeroMapUL edgeToCSRMap; /*!< \brief Map edges to CSR entries referenced by them (i,j) and (j,i). */

  /*--- Edge and element colorings. ---*/

  CCompressedSparsePatternUL edgeColoring, /*!< \brief Edge coloring structure for thread-based parallelization. */
      elemColoring;                        /*!< \brief Element coloring structure for thread-based parallelization. */
  unsigned long edgeColorGroupSize{1};     /*!< \brief Size of the edge groups within each color. */
  unsigned long elemColorGroupSize{1};     /*!< \brief Size of the element groups within each color. */

  ColMajorMatrix<uint8_t> CoarseGridColor_; /*!< \brief Coarse grid levels, colorized. */

 public:
  /*!< \brief Linelets (mesh lines perpendicular to stretching direction). */
  struct CLineletInfo {
    /*!< \brief Detect isotropic mesh region. */
    static passivedouble ALPHA_ISOTROPIC() { return 0.8; }
    enum : unsigned long { MAX_LINELET_POINTS = 32 }; /*!< \brief Maximum points per linelet. */

    std::vector<std::vector<unsigned long>> linelets; /*!< \brief Point indices for each linelet. */

    /*!< \brief Index of the linelet of each point ("linelets" transfered to points). */
    std::vector<unsigned> lineletIdx;

    /*!< \brief Signals that a point is not on a linelet. */
    enum : unsigned { NO_LINELET = std::numeric_limits<unsigned>::max() };

    /*!< \brief Coloring for OpenMP parallelization, "linelets" is sorted by color. */
    std::vector<unsigned long> colorOffsets;

    std::vector<uint8_t> lineletColor; /*!< \brief Coloring transfered to points, for visualization. */
  };

 protected:
  mutable CLineletInfo lineletInfo;

 public:
  /*--- Main geometric elements of the grid. ---*/

  CPrimalGrid** elem{nullptr};   /*!< \brief Element vector (primal grid information). */
  CPrimalGrid*** bound{nullptr}; /*!< \brief Boundary vector (primal grid information). */
  CPoint* nodes{nullptr};        /*!< \brief Node vector (dual grid information). */
  CEdge* edges{nullptr};         /*!< \brief Edge vector (dual grid information). */
  CVertex*** vertex{nullptr};    /*!< \brief Boundary Vertex vector (dual grid information). */
  CTurboVertex**** turbovertex{
      nullptr}; /*!< \brief Boundary Vertex vector ordered for turbomachinery calculation(dual grid information). */
  unsigned long* nVertex{nullptr};     /*!< \brief Number of vertex for each marker. */
  unsigned long* nElem_Bound{nullptr}; /*!< \brief Number of elements of the boundary. */
  string* Tag_to_Marker{nullptr};      /*!< \brief Names of boundary markers. */
  vector<bool>
      bound_is_straight; /*!< \brief Bool if boundary-marker is straight(2D)/plane(3D) for each local marker. */
  vector<su2double> SurfaceAreaCfgFile; /*!< \brief Total Surface area for all markers. */

  /*--- Partitioning-specific variables ---*/

  unordered_map<unsigned long, unsigned long>
      Global_to_Local_Elem;         /*!< \brief Mapping of global to local index for elements. */
  unsigned long* beg_node{nullptr}; /*!< \brief Array containing the first node on each rank due to a linear
                                       partitioning by global index. */
  unsigned long* end_node{
      nullptr}; /*!< \brief Array containing the last node on each rank due to a linear partitioning by global index. */
  unsigned long* nPointLinear{nullptr};     /*!< \brief Array containing the total number of nodes on each rank due to a
                                               linear partioning by global index. */
  unsigned long* nPointCumulative{nullptr}; /*!< \brief Cumulative storage array containing the total number of points
                                               on all prior ranks in the linear partitioning. */

  /*--- Data structures for point-to-point MPI communications. ---*/

  int maxCountPerPoint{0}; /*!< \brief Maximum number of pieces of data sent per vertex in point-to-point comms. */
  int nP2PSend{0};         /*!< \brief Number of sends during point-to-point comms. */
  int nP2PRecv{0};         /*!< \brief Number of receives during point-to-point comms. */
  int* nPoint_P2PSend{
      nullptr}; /*!< \brief Data structure holding number of vertices for each send in point-to-point comms. */
  int* nPoint_P2PRecv{
      nullptr}; /*!< \brief Data structure holding number of vertices for each recv in point-to-point comms. */
  int* Neighbors_P2PSend{
      nullptr}; /*!< \brief Data structure holding the ranks of the neighbors for point-to-point send comms. */
  int* Neighbors_P2PRecv{
      nullptr}; /*!< \brief Data structure holding the ranks of the neighbors for point-to-point recv comms. */
  map<int, int> P2PSend2Neighbor; /*!< \brief Data structure holding the reverse mapping of the ranks of the neighbors
                                     for point-to-point send comms. */
  map<int, int> P2PRecv2Neighbor; /*!< \brief Data structure holding the reverse mapping of the ranks of the neighbors
                                     for point-to-point recv comms. */
  unsigned long *Local_Point_P2PSend{nullptr}, /*!< \brief Data structure holding the local index of all vertices to be
                                                  sent in point-to-point comms. */
      *Local_Point_P2PRecv{nullptr}; /*!< \brief Data structure holding the local index of all vertices to be received
                                        in point-to-point comms. */
  su2double* bufD_P2PRecv{nullptr};  /*!< \brief Data structure for su2double point-to-point receive. */
  su2double* bufD_P2PSend{nullptr};  /*!< \brief Data structure for su2double point-to-point send. */
  unsigned short* bufS_P2PRecv{nullptr};  /*!< \brief Data structure for unsigned long point-to-point receive. */
  unsigned short* bufS_P2PSend{nullptr};  /*!< \brief Data structure for unsigned long point-to-point send. */
  SU2_MPI::Request* req_P2PSend{nullptr}; /*!< \brief Data structure for point-to-point send requests. */
  SU2_MPI::Request* req_P2PRecv{nullptr}; /*!< \brief Data structure for point-to-point recv requests. */

  /*--- Data structures for periodic communications. ---*/

  int maxCountPerPeriodicPoint{0}; /*!< \brief Maximum number of pieces of data sent per vertex in periodic comms. */
  int nPeriodicSend{0};            /*!< \brief Number of sends during periodic comms. */
  int nPeriodicRecv{0};            /*!< \brief Number of receives during periodic comms. */
  int* nPoint_PeriodicSend{
      nullptr}; /*!< \brief Data structure holding number of vertices for each send in periodic comms. */
  int* nPoint_PeriodicRecv{
      nullptr}; /*!< \brief Data structure holding number of vertices for each recv in periodic comms. */
  int* Neighbors_PeriodicSend{
      nullptr}; /*!< \brief Data structure holding the ranks of the neighbors for periodic send comms. */
  int* Neighbors_PeriodicRecv{
      nullptr}; /*!< \brief Data structure holding the ranks of the neighbors for periodic recv comms. */
  map<int, int> PeriodicSend2Neighbor; /*!< \brief Data structure holding the reverse mapping of the ranks of the
                                          neighbors for periodic send comms. */
  map<int, int> PeriodicRecv2Neighbor; /*!< \brief Data structure holding the reverse mapping of the ranks of the
                                          neighbors for periodic recv comms. */
  unsigned long *Local_Point_PeriodicSend{
      nullptr}, /*!< \brief Data structure holding the local index of all vertices to be sent in periodic comms. */
      *Local_Point_PeriodicRecv{nullptr},  /*!< \brief Data structure holding the local index of all vertices to be
                                              received in periodic comms. */
      *Local_Marker_PeriodicSend{nullptr}, /*!< \brief Data structure holding the local index of the periodic marker for
                                              a particular vertex to be sent in periodic comms. */
      *Local_Marker_PeriodicRecv{nullptr}; /*!< \brief Data structure holding the local index of the periodic marker for
                                              a particular vertex to be received in periodic comms. */
  su2double* bufD_PeriodicRecv{nullptr};   /*!< \brief Data structure for su2double periodic receive. */
  su2double* bufD_PeriodicSend{nullptr};   /*!< \brief Data structure for su2double periodic send. */
  unsigned short* bufS_PeriodicRecv{nullptr};  /*!< \brief Data structure for unsigned long periodic receive. */
  unsigned short* bufS_PeriodicSend{nullptr};  /*!< \brief Data structure for unsigned long periodic send. */
  SU2_MPI::Request* req_PeriodicSend{nullptr}; /*!< \brief Data structure for periodic send requests. */
  SU2_MPI::Request* req_PeriodicRecv{nullptr}; /*!< \brief Data structure for periodic recv requests. */

  /*--- Mesh quality metrics. ---*/

  vector<su2double> Orthogonality; /*!< \brief Measure of dual CV orthogonality angle (0 to 90 deg., 90 being best). */
  vector<su2double> Aspect_Ratio;  /*!< \brief Measure of dual CV aspect ratio (max face area / min face area).  */
  vector<su2double>
      Volume_Ratio; /*!< \brief Measure of dual CV volume ratio (max sub-element volume / min sub-element volume). */

  const ColMajorMatrix<uint8_t>& CoarseGridColor = CoarseGridColor_; /*!< \brief Coarse grid levels, colorized. */

  /*!
   * \brief Constructor of the class.
   */
  CGeometry();

  /*!
   * \brief Constructor of the class.
   */
  CGeometry(CConfig* config, unsigned short nDim);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CGeometry();

  /*!
   * \brief Routine to launch non-blocking recvs only for all periodic communications.
   * \brief Routine to set up persistent data structures for point-to-point MPI communications.
   * \note This routine is called by any class that has loaded data into the generic communication buffers.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void PreprocessP2PComms(CGeometry* geometry, CConfig* config);

  /*!
   * \brief Routine to allocate buffers for point-to-point MPI communications. Also called to dynamically reallocate if
   * not enough memory is found for comms during runtime. \param[in] val_countPerPoint - Maximum count of the data type
   * per vertex in point-to-point comms, e.g., nPrimvarGrad*nDim.
   */
  void AllocateP2PComms(unsigned short val_countPerPoint);

  /*!
   * \brief Routine to launch non-blocking recvs only for all point-to-point communication with neighboring partitions.
   * \note This routine is called by any class that has loaded data into the generic communication buffers.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be communicated.
   * \param[in] countPerPoint - Number of variables per point.
   * \param[in] val_reverse - Boolean controlling forward or reverse communication between neighbors.
   */
  void PostP2PRecvs(CGeometry* geometry, const CConfig* config, unsigned short commType, unsigned short countPerPoint,
                    bool val_reverse) const;

  /*!
   * \brief Routine to launch a single non-blocking send once the buffer is loaded for a point-to-point commucation.
   * \note This routine is called by any class that has loaded data into the generic communication buffers.
   * \param[in] geometry     - Geometrical definition of the problem.
   * \param[in] config       - Definition of the particular problem.
   * \param[in] commType     - Enumerated type for the quantity to be communicated.
   * \param[in] countPerPoint - Number of variables per point.
   * \param[in] val_iMessage - Index of the message in the order they are stored.
   * \param[in] val_reverse  - Boolean controlling forward or reverse communication between neighbors.
   */
  void PostP2PSends(CGeometry* geometry, const CConfig* config, unsigned short commType, unsigned short countPerPoint,
                    int val_iMessage, bool val_reverse) const;

  /*!
   * \brief Routine to set up persistent data structures for periodic communications.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void PreprocessPeriodicComms(CGeometry* geometry, CConfig* config);

  /*!
   * \brief Routine to allocate buffers for periodic communications. Also called to dynamically reallocate if not enough
   * memory is found for comms during runtime. \param[in] val_countPerPeriodicPoint - Maximum count of the data type per
   * vertex in periodic comms, e.g., nPrimvarGrad*nDim.
   */
  void AllocatePeriodicComms(unsigned short val_countPerPeriodicPoint);

  /*!
   * \brief Routine to launch non-blocking recvs only for all periodic communication with neighboring partitions.
   * \note This routine is called by any class that has loaded data into the generic communication buffers.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be communicated.
   * \param[in] countPerPeriodicPoint - Number of variables per point.
   */
  void PostPeriodicRecvs(CGeometry* geometry, const CConfig* config, unsigned short commType,
                         unsigned short countPerPeriodicPoint);

  /*!
   * \brief Routine to launch a single non-blocking send once the buffer is loaded for a periodic commucation.
   * \note This routine is called by any class that has loaded data into the generic communication buffers.
   * \param[in] geometry     - Geometrical definition of the problem.
   * \param[in] config       - Definition of the particular problem.
   * \param[in] commType     - Enumerated type for the quantity to be communicated.
   * \param[in] countPerPeriodicPoint - Number of variables per point.
   * \param[in] val_iMessage - Index of the message in the order they are stored.
   */
  void PostPeriodicSends(CGeometry* geometry, const CConfig* config, unsigned short commType,
                         unsigned short countPerPeriodicPoint, int val_iMessage) const;

  /*!
   * \brief Helper function to define the type and number of variables per point for each communication type.
   * \param[in] config - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be communicated.
   * \param[out] COUNT_PER_POINT - Number of communicated variables per point.
   * \param[out] MPI_TYPE - Enumerated type for the datatype of the quantity to be communicated.
   */
  void GetCommCountAndType(const CConfig* config, unsigned short commType, unsigned short& COUNT_PER_POINT,
                           unsigned short& MPI_TYPE) const;

  /*!
   * \brief Routine to load a geometric quantity into the data structures for MPI point-to-point communication and to
   *        launch non-blocking sends and recvs for all point-to-point communication with neighboring partitions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be communicated.
   */
  void InitiateComms(CGeometry* geometry, const CConfig* config, unsigned short commType) const;

  /*!
   * \brief Routine to complete the set of non-blocking communications launched by InitiateComms() and unpacking of the
   * data into the geometry class. \param[in] geometry - Geometrical definition of the problem. \param[in] config   -
   * Definition of the particular problem. \param[in] commType - Enumerated type for the quantity to be unpacked.
   */
  void CompleteComms(CGeometry* geometry, const CConfig* config, unsigned short commType);

  /*!
   * \brief Get number of coordinates.
   * \return Number of coordinates.
   */
  inline unsigned short GetnDim() const { return nDim; }

  /*!
   * \brief Get number of zones.
   * \return Number of zones.
   */
  inline unsigned short GetnZone() const { return nZone; }

  /*!
   * \brief Get number of points.
   * \return Number of points.
   */
  inline unsigned long GetnPoint() const { return nPoint; }

  /*!
   * \brief Get number of real points (that belong to the domain).
   * \return Number of real points.
   */
  inline unsigned long GetnPointDomain() const { return nPointDomain; }

  /*!
   * \brief Retrieve total number of nodes in a simulation across all processors (including halos).
   * \return Total number of nodes in a simulation across all processors (including halos).
   */
  inline unsigned long GetGlobal_nPoint() const { return Global_nPoint; }

  /*!
   * \brief Retrieve total number of nodes in a simulation across all processors (excluding halos).
   * \return Total number of nodes in a simulation across all processors (excluding halos).
   */
  inline unsigned long GetGlobal_nPointDomain() const { return Global_nPointDomain; }

  /*!
   * \brief Get number of elements.
   * \return Number of elements.
   */
  inline unsigned long GetnElem() const { return nElem; }

  /*!
   * \brief Get number of edges.
   * \return Number of edges.
   */
  inline unsigned long GetnEdge() const { return nEdge; }

  /*!
   * \brief Get number of markers.
   * \return Number of markers.
   */
  inline unsigned short GetnMarker() const { return nMarker; }

  /*!
   * \brief Get number of vertices.
   * \param[in] val_marker - Marker of the boundary.
   * \return Number of vertices.
   */
  inline const su2double* GetSpanWiseValue(unsigned short val_marker) const { return SpanWiseValue[val_marker - 1]; }

  /*!
   * \brief Get number of vertices.
   * \param[in] val_marker - Marker of the boundary.
   * \return Number of vertices.
   */
  inline unsigned long GetnVertex(unsigned short val_marker) const { return nVertex[val_marker]; }

  /*!
   * \brief Get number of span wise section.
   * \param[in] marker_flag - flag of the turbomachinery boundary.
   * \return Number of span wise section.
   */
  inline unsigned short GetnSpanWiseSections(unsigned short marker_flag) const {
    return nSpanWiseSections[marker_flag - 1];
  }

  /*!
   * \brief Get number of vertices.
   * \param[in] val_marker - Marker of the boundary.
   * \return Number of vertices.
   */
  inline unsigned long GetnVertexSpan(unsigned short val_marker, unsigned short val_span) const {
    return nVertexSpan[val_marker][val_span];
  }

  /*!
   * \brief Get number of frequencies per span for NRBC.
   * \param[in] val_marker - Marker of the boundary.
   * \return Number of frequencies for NRBC.
   */
  inline unsigned long GetnFreqSpan(unsigned short val_marker, unsigned short val_span) const {
    return (nTotVertexSpan[val_marker][val_span] / 2 - 1);
  }

  /*!
   * \brief Get number of vertices.
   * \param[in] val_marker - Marker of the boundary.
   * \return Number of vertices.
   */
  inline unsigned long GetnVertexSpanMax(unsigned short marker_flag) const { return nVertexSpanMax[marker_flag]; }

  /*!
   * \brief Get number of max frequencies for initializing the Fourier Coefficient for NR BC.
   * \param[in] marker_flag - Marker of the boundary.
   * \return Number of frequencies.
   */
  inline unsigned long GetnFreqSpanMax(unsigned short marker_flag) const {
    return (nVertexSpanMax[marker_flag] / 2 - 1);
  }

  /*!
   * \brief Get number of vertices.
   * \param[in] val_marker - Marker of the boundary.
   * \return Number of vertices.
   */
  inline void SetnVertexSpanMax(unsigned short marker_flag, unsigned long nVertMax) {
    nVertexSpanMax[marker_flag] = nVertMax;
  }

  /*!
   * \brief Get the edge index from using the nodes of the edge.
   * \param[in] first_point - First point of the edge.
   * \param[in] second_point - Second point of the edge.
   * \param[in] error - Throw error if edge does not exist.
   * \return Index of the edge.
   */
  inline long FindEdge(unsigned long first_point, unsigned long second_point, bool error = true) const {
    for (unsigned short iNode = 0; iNode < nodes->GetnPoint(first_point); iNode++) {
      auto iPoint = nodes->GetPoint(first_point, iNode);
      if (iPoint == second_point) return nodes->GetEdge(first_point, iNode);
    }
    if (error) {
      char buf[100];
      SPRINTF(buf, "Can't find the edge that connects %lu and %lu.", first_point, second_point);
      SU2_MPI::Error(buf, CURRENT_FUNCTION);
    }
    return -1;
  }

  /*!
   * \brief Get the edge index from using the nodes of the edge.
   * \param[in] first_point - First point of the edge.
   * \param[in] second_point - Second point of the edge.
   * \return Index of the edge.
   */
  inline bool CheckEdge(unsigned long first_point, unsigned long second_point) const {
    return FindEdge(first_point, second_point, false) >= 0;
  }

  /*!
   * \brief Create a file for testing the geometry.
   */
  void TestGeometry() const;

  /*!
   * \brief A virtual member.
   * \param[in] val_nmarker - Number of markers.
   */
  inline void SetnMarker(unsigned short val_nmarker) { nMarker = val_nmarker; }

  /*!
   * \brief Set the number of dimensions of the problem.
   * \param[in] val_nDim - Number of dimensions.
   */
  inline void SetnDim(unsigned short val_nDim) { nDim = val_nDim; }

  /*!
   * \brief Get the index of a marker.
   * \param[in] val_marker - Marker of the boundary.
   * \return Index of the marker in the grid defintion.
   */
  inline string GetMarker_Tag(unsigned short val_marker) const { return Tag_to_Marker[val_marker]; }

  /*!
   * \brief Set index of a marker.
   * \param[in] val_marker - Marker of the boundary.
   * \param[in] val_index - Index of the marker.
   */
  inline void SetMarker_Tag(unsigned short val_marker, string val_index) {
    Tag_to_Marker[val_marker] = std::move(val_index);
  }

  /*!
   * \brief Set the number of boundary elements.
   * \param[in] val_marker - Marker of the boundary.
   * \param[in] val_nelem_bound - Number of boundary elements.
   */
  inline void SetnElem_Bound(unsigned short val_marker, unsigned long val_nelem_bound) {
    nElem_Bound[val_marker] = val_nelem_bound;
  }

  /*!
   * \brief Set the number of grid points.
   * \param[in] val_npoint - Number of grid points.
   */
  inline void SetnPoint(unsigned long val_npoint) { nPoint = val_npoint; }

  /*!
   * \brief Set the number of grid points in the domain.
   * \param[in] val_npoint - Number of grid points in the domain.
   */
  inline void SetnPointDomain(unsigned long val_npoint) { nPointDomain = val_npoint; }

  /*!
   * \brief Set the value of the total number of points globally in the simulation.
   * \param[in] val_global_npoint - Global number of points in the mesh (excluding halos).
   */
  void SetGlobal_nPointDomain(unsigned long val_global_npoint) { Global_nPointDomain = val_global_npoint; }

  /*!
   * \brief Set the number of grid elements.
   * \param[in] val_nelem - Number of grid elements.
   */
  inline void SetnElem(unsigned long val_nelem) { nElem = val_nelem; }

  /*!
   * \brief Get the number of boundary elements.
   * \param[in] val_marker - Marker of the boundary.
   */
  inline unsigned long GetnElem_Bound(unsigned short val_marker) const { return nElem_Bound[val_marker]; }

  /*!
   * \brief Get the number of elements in vtk fortmat.
   */
  inline unsigned long GetMax_GlobalPoint() const { return Max_GlobalPoint; }

  /*!
   * \brief Finds face of element.
   * \param[in] first_elem - Identification of the first element.
   * \param[in] second_elem - Identification of the second element.
   * \param[in] face_first_elem - Index of the common face for the first element.
   * \param[in] face_second_elem - Index of the common face for the second element.
   */
  inline virtual bool FindFace(unsigned long first_elem, unsigned long second_elem, unsigned short& face_first_elem,
                               unsigned short& face_second_elem) {
    return false;
  }

  /*!
   * \brief Sets area to be positive in Z direction.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetPositive_ZArea(CConfig* config) {}

  /*!
   * \brief Set connectivity between points.
   */
  inline virtual void SetPoint_Connectivity() {}

  /*!
   * \brief Orders the RCM.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetRCM_Ordering(CConfig* config) {}

  /*!
   * \brief Connects elements  .
   */
  inline virtual void SetElement_Connectivity() {}

  /*!
   * \brief Sets the edges of an elemment.
   */
  void SetEdges();

  /*!
   * \brief Sets the faces of an element..
   */
  void SetFaces();

  /*!
   * \brief Sets the boundary volume.
   */
  inline virtual void SetBoundVolume() {}

  /*!
   * \brief Sets the vertices.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetVertex(const CConfig* config) {}

  /*!
   * \brief Computes the N span.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iZone - Zone of the problem
   * \param[in] marker_flag - Marker being used
   * \param[in] allocate
   */
  inline virtual void ComputeNSpan(CConfig* config, unsigned short val_iZone, unsigned short marker_flag,
                                   bool allocate) {}

  /*!
   * \brief Set vertices for turbomachinery problems.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iZone - Zone of the problem
   * \param[in] marker_flag - Marker being used
   * \param[in] allocate
   */
  inline virtual void SetTurboVertex(CConfig* config, unsigned short val_iZone, unsigned short marker_flag,
                                     bool allocate) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iZone - Zone of the problem
   * \param[in] marker_flag - Marker being used
   */
  inline virtual void UpdateTurboVertex(CConfig* config, unsigned short val_iZone, unsigned short marker_flag) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iZone - Zone of the problem
   * \param[in] marker_flag - Marker being used
   * \param[in] allocate
   */
  inline virtual void SetAvgTurboValue(CConfig* config, unsigned short val_iZone, unsigned short marker_flag,
                                       bool allocate) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   * \param[in] allocate
   */
  inline virtual void GatherInOutAverageValues(CConfig* config, bool allocate) {}

  /*!
   * \brief Set max length.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetMaxLength(CConfig* config) {}

  /*!
   * \brief  Sets control volume.
   * \param[in] config - Definition of the particular problem.
   * \param[in] action - Allocate or not the new elements.
   */
  inline virtual void SetControlVolume(CConfig* config, unsigned short action) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void VisualizeControlVolume(const CConfig* config) const {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void MatchActuator_Disk(const CConfig* config) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void MatchPeriodic(const CConfig* config, unsigned short val_periodic) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   * \param[in] action - Allocate or not the new elements.
   */
  inline virtual void SetBoundControlVolume(const CConfig* config, unsigned short action) {}

  /*!
   * \brief A virtual member.
   * \param[in] config_filename - Name of the file where the tecplot information is going to be stored.
   */
  inline virtual void SetTecPlot(char config_filename[MAX_STRING_SIZE], bool new_file) {}

  /*!
   * \brief A virtual member.
   * \param[in] mesh_filename - Name of the file where the tecplot information is going to be stored.
   * \param[in] new_file - Boolean to decide if aopen a new file or add to a old one
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetBoundTecPlot(char mesh_filename[MAX_STRING_SIZE], bool new_file, CConfig* config) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Check_IntElem_Orientation(const CConfig* config) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Check_BoundElem_Orientation(const CConfig* config) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetColorGrid(CConfig* config) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetColorGrid_Parallel(const CConfig* config) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetColorFEMGrid_Parallel(CConfig* config) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void DivideConnectivity(CConfig* config, unsigned short Elem_Type) {}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_domain - Number of domains for parallelization purposes.
   */
  inline virtual void SetSendReceive(const CConfig* config) {}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_domain - Number of domains for parallelization purposes.
   */
  inline virtual void SetBoundaries(CConfig* config) {}

  /*!
   * \brief A virtual member.
   * \param[in] fine_grid - Geometrical definition of the problem.
   */
  inline virtual void SetCoord(const CGeometry* fine_grid) {}

  /*!
   * \brief A virtual member.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] val_marker - Index of the boundary marker.
   */
  inline virtual void SetMultiGridWallHeatFlux(const CGeometry* fine_grid, unsigned short val_marker) {}

  /*!
   * \brief A virtual member.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] val_marker - Index of the boundary marker.
   */
  inline virtual void SetMultiGridWallTemperature(const CGeometry* fine_grid, unsigned short val_marker) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_nSmooth - Number of smoothing iterations.
   * \param[in] val_smooth_coeff - Relaxation factor.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetCoord_Smoothing(unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig* config) {}

  /*!
   * \brief A virtual member.
   * \param[in] fine_grid - Geometrical definition of the child grid (for multigrid).
   */
  inline virtual void SetPoint_Connectivity(const CGeometry* fine_grid) {}

  /*!
   * \brief A virtual member.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetVertex(const CGeometry* fine_grid, const CConfig* config) {}

  /*!
   * \brief A virtual member.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] action - Allocate or not the new elements.
   */
  inline virtual void SetControlVolume(const CGeometry* fine_grid, unsigned short action) {}

  /*!
   * \brief A virtual member.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] action - Allocate or not the new elements.
   */
  inline virtual void SetBoundControlVolume(const CGeometry* fine_grid, unsigned short action) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetBoundSensitivity(CConfig* config) {}

  /*!
   * \brief Set the data containers for customized boundary conditions.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCustomBoundary(CConfig* config);

  /*!
   * \brief Set cartesian grid velocity based on rotational speed and axis.
   * \param[in] config - Definition of the particular problem.
   * \param[in] print - Display information on screen.
   */
  void SetRotationalVelocity(const CConfig* config, bool print = false);

  /*!
   * \brief Set the rotational velocity of the points on the shroud markers to 0.
   * \param[in] config - Definition of the particular problem.
   */
  void SetShroudVelocity(const CConfig* config);

  /*!
   * \brief Set the translational velocity at each node.
   * \param[in] config - Definition of the particular problem.
   * \param[in] print - Display information on screen.
   */
  void SetTranslationalVelocity(const CConfig* config, bool print = false);

  /*!
   * \brief Set the translational/rotational velocity for all moving walls.
   * \param[in] config - Definition of the particular problem.
   * \param[in] print - Display information on screen.
   */
  void SetWallVelocity(const CConfig* config, bool print = false);

  /*!
   * \brief Set the grid velocity via finite differencing at each node.
   * \param[in] config - Definition of the particular problem.
   */
  void SetGridVelocity(const CConfig* config);

  /*!
   * \brief A virtual member.
   * \param[in] fine_grid - Geometry of the fine mesh.
   */
  inline virtual void SetRestricted_GridVelocity(const CGeometry* fine_grid) {}

  /*!
   * \brief Compute the surface area of all global markers.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeSurfaceAreaCfgFile(const CConfig* config);

  /*!
   * \brief Get global Surface Area to a local marker.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Local surface marker.
   * \return Global Surface Area to the local marker
   */
  su2double GetSurfaceArea(const CConfig* config, unsigned short val_marker) const;

  /*!
   * \brief Check if a boundary is straight(2D) / plane(3D) for EULER_WALL and SYMMETRY_PLANE
   *        only and store the information in bound_is_straight. For all other boundary types
   *        this will return false and could therfore be wrong. Used ultimately for BC_Slip_Wall.
   * \param[in] config - Definition of the particular problem.
   * \param[in] print_on_screen - Boolean whether to print result on screen.
   */
  void ComputeSurf_Straightness(CConfig* config, bool print_on_screen);

  /*!
   * \brief Find and store all vertices on a sharp corner in the geometry.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeSurf_Curvature(CConfig* config);

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeAirfoil_Section(su2double* Plane_P0, su2double* Plane_Normal, su2double MinXCoord, su2double MaxXCoord,
                              su2double MinYCoord, su2double MaxYCoord, su2double MinZCoord, su2double MaxZCoord,
                              const su2double* FlowVariable, vector<su2double>& Xcoord_Airfoil,
                              vector<su2double>& Ycoord_Airfoil, vector<su2double>& Zcoord_Airfoil,
                              vector<su2double>& Variable_Airfoil, bool original_surface, CConfig* config);

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_MaxThickness(su2double* Plane_P0, su2double* Plane_Normal, CConfig* config,
                                         vector<su2double>& Xcoord_Airfoil, vector<su2double>& Ycoord_Airfoil,
                                         vector<su2double>& Zcoord_Airfoil) {
    return 0.0;
  }

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_Twist(su2double* Plane_P0, su2double* Plane_Normal, vector<su2double>& Xcoord_Airfoil,
                                  vector<su2double>& Ycoord_Airfoil, vector<su2double>& Zcoord_Airfoil) {
    return 0.0;
  }

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_Chord(su2double* Plane_P0, su2double* Plane_Normal, vector<su2double>& Xcoord_Airfoil,
                                  vector<su2double>& Ycoord_Airfoil, vector<su2double>& Zcoord_Airfoil) {
    return 0.0;
  }

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_Width(su2double* Plane_P0, su2double* Plane_Normal, vector<su2double>& Xcoord_Airfoil,
                                  vector<su2double>& Ycoord_Airfoil, vector<su2double>& Zcoord_Airfoil) {
    return 0.0;
  }

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_WaterLineWidth(su2double* Plane_P0, su2double* Plane_Normal, CConfig* config,
                                           vector<su2double>& Xcoord_Airfoil, vector<su2double>& Ycoord_Airfoil,
                                           vector<su2double>& Zcoord_Airfoil) {
    return 0.0;
  }

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_Height(su2double* Plane_P0, su2double* Plane_Normal, vector<su2double>& Xcoord_Airfoil,
                                   vector<su2double>& Ycoord_Airfoil, vector<su2double>& Zcoord_Airfoil) {
    return 0.0;
  }

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_LERadius(su2double* Plane_P0, su2double* Plane_Normal, vector<su2double>& Xcoord_Airfoil,
                                     vector<su2double>& Ycoord_Airfoil, vector<su2double>& Zcoord_Airfoil) {
    return 0.0;
  }

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_Thickness(su2double* Plane_P0, su2double* Plane_Normal, su2double Location, CConfig* config,
                                      vector<su2double>& Xcoord_Airfoil, vector<su2double>& Ycoord_Airfoil,
                                      vector<su2double>& Zcoord_Airfoil, su2double& ZLoc) {
    return 0.0;
  }

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_Area(su2double* Plane_P0, su2double* Plane_Normal, CConfig* config,
                                 vector<su2double>& Xcoord_Airfoil, vector<su2double>& Ycoord_Airfoil,
                                 vector<su2double>& Zcoord_Airfoil) {
    return 0.0;
  }

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_Length(su2double* Plane_P0, su2double* Plane_Normal, CConfig* config,
                                   vector<su2double>& Xcoord_Airfoil, vector<su2double>& Ycoord_Airfoil,
                                   vector<su2double>& Zcoord_Airfoil) {
    return 0.0;
  }

  /*!
   * \brief A virtual member.
   */
  virtual void Compute_Wing_LeadingTrailing(su2double* LeadingEdge, su2double* TrailingEdge, su2double* Plane_P0,
                                            su2double* Plane_Normal, vector<su2double>& Xcoord_Airfoil,
                                            vector<su2double>& Ycoord_Airfoil, vector<su2double>& Zcoord_Airfoil) {}

  /*!
   * \brief A virtual member.
   */
  virtual void Compute_Fuselage_LeadingTrailing(su2double* LeadingEdge, su2double* TrailingEdge, su2double* Plane_P0,
                                                su2double* Plane_Normal, vector<su2double>& Xcoord_Airfoil,
                                                vector<su2double>& Ycoord_Airfoil, vector<su2double>& Zcoord_Airfoil) {}

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_Dihedral(su2double* LeadingEdge_im1, su2double* TrailingEdge_im1, su2double* LeadingEdge_i,
                                     su2double* TrailingEdge_i) {
    return 0.0;
  }

  /*!
   * \brief A virtual member.
   */
  virtual su2double Compute_Curvature(su2double* LeadingEdge_im1, su2double* TrailingEdge_im1, su2double* LeadingEdge_i,
                                      su2double* TrailingEdge_i, su2double* LeadingEdge_ip1,
                                      su2double* TrailingEdge_ip1) {
    return 0.0;
  }

  /*!
   * \brief A virtual member.
   */
  virtual void Compute_Wing(CConfig* config, bool original_surface, su2double& Wing_Volume,
                            su2double& Wing_MinMaxThickness, su2double& Wing_MaxMaxThickness, su2double& Wing_MinChord,
                            su2double& Wing_MaxChord, su2double& Wing_MinLERadius, su2double& Wing_MaxLERadius,
                            su2double& Wing_MinToC, su2double& Wing_MaxToC, su2double& Wing_ObjFun_MinToC,
                            su2double& Wing_MaxTwist, su2double& Wing_MaxCurvature, su2double& Wing_MaxDihedral) {}

  /*!
   * \brief A virtual member.
   */
  virtual void Compute_Fuselage(CConfig* config, bool original_surface, su2double& Fuselage_Volume,
                                su2double& Fuselage_WettedArea, su2double& Fuselage_MinWidth,
                                su2double& Fuselage_MaxWidth, su2double& Fuselage_MinWaterLineWidth,
                                su2double& Fuselage_MaxWaterLineWidth, su2double& Fuselage_MinHeight,
                                su2double& Fuselage_MaxHeight, su2double& Fuselage_MaxCurvature) {}

  /*!
   * \brief A virtual member.
   */
  virtual void Compute_Nacelle(CConfig* config, bool original_surface, su2double& Nacelle_Volume,
                               su2double& Nacelle_MinMaxThickness, su2double& Nacelle_MaxMaxThickness,
                               su2double& Nacelle_MinChord, su2double& Nacelle_MaxChord, su2double& Nacelle_MinLERadius,
                               su2double& Nacelle_MaxLERadius, su2double& Nacelle_MinToC, su2double& Nacelle_MaxToC,
                               su2double& Nacelle_ObjFun_MinToC, su2double& Nacelle_MaxTwist) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void FindNormal_Neighbor(const CConfig* config) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetGlobal_to_Local_Point() {}

  /*!
   * \brief A virtual member.
   * \param[in] val_ipoint - Global point.
   * \return Local index that correspond with the global index.
   */
  inline virtual long GetGlobal_to_Local_Point(unsigned long val_ipoint) const { return 0; }

  /*!
   * \brief Retrieve total number of elements in a simulation across all processors.
   * \return Total number of elements in a simulation across all processors.
   */
  inline unsigned long GetGlobal_nElem() const { return Global_nElem; }

  /*!
   * \brief  Retrieve total number of elements in a simulation across all processors (excluding halos).
   * \return Total number of elements in a simulation across all processors (excluding halos).
   */
  inline unsigned long GetGlobal_nElemDomain() const { return Global_nElemDomain; }

  /*!
   * \brief Retrieve total number of triangular elements in a simulation across all processors.
   * \return Total number of line elements in a simulation across all processors.
   */
  inline unsigned long GetGlobal_nElemLine() const { return Global_nelem_edge; }

  /*!
   * \brief Retrieve total number of triangular elements in a simulation across all processors.
   * \return Total number of triangular elements in a simulation across all processors.
   */
  inline unsigned long GetGlobal_nElemTria() const { return Global_nelem_triangle; }

  /*!
   * \brief Retrieve total number of quadrilateral elements in a simulation across all processors.
   * \return Total number of quadrilateral elements in a simulation across all processors.
   */
  inline unsigned long GetGlobal_nElemQuad() const { return Global_nelem_quad; }

  /*!
   * \brief Retrieve total number of tetrahedral elements in a simulation across all processors.
   * \return Total number of tetrahedral elements in a simulation across all processors.
   */
  inline unsigned long GetGlobal_nElemTetr() const { return Global_nelem_tetra; }

  /*!
   * \brief Retrieve total number of hexahedral elements in a simulation across all processors.
   * \return Total number of hexahedral elements in a simulation across all processors.
   */
  inline unsigned long GetGlobal_nElemHexa() const { return Global_nelem_hexa; }

  /*!
   * \brief Retrieve total number of prism elements in a simulation across all processors.
   * \return Total number of prism elements in a simulation across all processors.
   */
  inline unsigned long GetGlobal_nElemPris() const { return Global_nelem_prism; }

  /*!
   * \brief Retrieve total number of pyramid elements in a simulation across all processors.
   * \return Total number of pyramid elements in a simulation across all processors.
   */
  inline unsigned long GetGlobal_nElemPyra() const { return Global_nelem_pyramid; }

  /*!
   * \brief Get number of triangular elements.
   * \return Number of line elements.
   */
  inline unsigned long GetnElemLine() const { return nelem_edge; }

  /*!
   * \brief Get number of triangular elements.
   * \return Number of triangular elements.
   */
  inline unsigned long GetnElemTria() const { return nelem_triangle; }

  /*!
   * \brief Get number of quadrilateral elements.
   * \return Number of quadrilateral elements.
   */
  inline unsigned long GetnElemQuad() const { return nelem_quad; }

  /*!
   * \brief Get number of tetrahedral elements.
   * \return Number of tetrahedral elements.
   */
  inline unsigned long GetnElemTetr() const { return nelem_tetra; }

  /*!
   * \brief Get number of hexahedral elements.
   * \return Number of hexahedral elements.
   */
  inline unsigned long GetnElemHexa() const { return nelem_hexa; }

  /*!
   * \brief Get number of prism elements.
   * \return Number of prism elements.
   */
  inline unsigned long GetnElemPris() const { return nelem_prism; }

  /*!
   * \brief Get number of pyramid elements.
   * \return Number of pyramid elements.
   */
  inline unsigned long GetnElemPyra() const { return nelem_pyramid; }

  /*!
   * \brief Get x coords of geometrical planes in the mesh
   */
  inline vector<vector<su2double>> GetXCoord() const { return Xcoord_plane; }

  /*!
   * \brief Get y coords of geometrical planes in the mesh
   */
  inline vector<vector<su2double>> GetYCoord() const { return Ycoord_plane; }

  /*!
   * \brief Get z coords of geometrical planes in the mesh
   */
  inline vector<vector<su2double>> GetZCoord() const { return Zcoord_plane; }

  /*!
   * \brief Get all points on a geometrical plane in the mesh
   */
  inline vector<vector<unsigned long>> GetPlanarPoints() const { return Plane_points; }

  /*!
   * \brief Compute the intersection between a segment and a plane.
   * \param[in] Segment_P0 - Definition of the particular problem.
   * \param[in] Segment_P1 - Definition of the particular problem.
   * \param[in] Plane_P0 - Definition of the particular problem.
   * \param[in] Plane_Normal - Definition of the particular problem.
   * \param[in] Intersection - Definition of the particular problem.
   * \return If the intersection has has been successful.
   */
  bool SegmentIntersectsPlane(const su2double* Segment_P0, const su2double* Segment_P1, su2double Variable_P0,
                              su2double Variable_P1, const su2double* Plane_P0, const su2double* Plane_Normal,
                              su2double* Intersection, su2double& Variable_Interp);

  /*!
   * \brief Ray Intersects Triangle (Moller and Trumbore algorithm)
   */
  bool RayIntersectsTriangle(const su2double orig[3], const su2double dir[3], const su2double vert0[3],
                             const su2double vert1[3], const su2double vert2[3], su2double* intersect);

  /*!
   * \brief Segment Intersects Triangle
   */
  bool SegmentIntersectsTriangle(su2double point0[3], const su2double point1[3], su2double vert0[3], su2double vert1[3],
                                 su2double vert2[3]);

  /*!
   * \brief Segment Intersects Line (for 2D FFD Intersection)
   */
  bool SegmentIntersectsLine(const su2double point0[2], const su2double point1[2], const su2double vert0[2],
                             const su2double vert1[2]);

  /*!
   * \brief Register the coordinates of the mesh nodes.
   */
  void RegisterCoordinates() const;

  /*!
   * \brief Update the multi-grid structure and the wall-distance.
   * \param geometry_container - Geometrical definition.
   * \param config - Config
   */
  static void UpdateGeometry(CGeometry** geometry_container, CConfig* config);

  /*!
   * \brief Update the multi-grid structure for the customized boundary conditions
   * \param geometry_container - Geometrical definition.
   * \param config - Definition of the particular problem.
   */
  void UpdateCustomBoundaryConditions(CGeometry** geometry_container, CConfig* config);

  /*!
   * \brief A virtual member.
   * \param config - Config
   */
  inline virtual void SetSensitivity(CConfig* config) {}

  /*!
   * \brief A virtual member.
   * \param config - Config
   */
  inline virtual void ReadUnorderedSensitivity(CConfig* config) {}

  /*!
   * \brief A virtual member.
   * \param iPoint - Point
   * \param iDim - Dimension
   */
  inline virtual su2double GetSensitivity(unsigned long iPoint, unsigned short iDim) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param iPoint - Point
   * \param iDim - Dimension
   * \param val - Value of the sensitivity
   */
  inline virtual void SetSensitivity(unsigned long iPoint, unsigned short iDim, su2double val) {}

  /*!
   * \brief Get the average normal at a specific span for a given marker in the turbomachinery reference of frame.
   * \param[in] val_marker - marker value.
   * \param[in] val_span - span value.
   * \return The span-wise averaged turbo normal.
   */
  inline const su2double* GetAverageTurboNormal(unsigned short val_marker, unsigned short val_span) const {
    return AverageTurboNormal[val_marker][val_span];
  }

  /*!
   * \brief Get the average normal at a specific span for a given marker.
   * \param[in] val_marker - marker value.
   * \param[in] val_span - span value.
   * \return The span-wise averaged normal.
   */
  inline const su2double* GetAverageNormal(unsigned short val_marker, unsigned short val_span) const {
    return AverageNormal[val_marker][val_span];
  }

  /*!
   * \brief Get the value of the total area for each span.
   * \param[in] val_marker - marker value.
   * \param[in] val_span - span value.
   * \return The span-wise area.
   */
  inline su2double GetSpanArea(unsigned short val_marker, unsigned short val_span) const {
    return SpanArea[val_marker][val_span];
  }

  /*!
   * \brief Get the value of the total area for each span.
   * \param[in] val_marker - marker value.
   * \param[in] val_span - span value.
   * \return The span-wise averaged turbo normal.
   */
  inline su2double GetTurboRadius(unsigned short val_marker, unsigned short val_span) const {
    return TurboRadius[val_marker][val_span];
  }

  /*!
   * \brief Get the value of the average tangential rotational velocity for each span.
   * \param[in] val_marker - marker value.
   * \param[in] val_span - span value.
   * \return The span-wise averaged tangential velocity.
   */
  inline su2double GetAverageTangGridVel(unsigned short val_marker, unsigned short val_span) const {
    return AverageTangGridVel[val_marker][val_span];
  }

  /*!
   * \brief Get the value of the inflow tangential velocity at each span.
   * \param[in] val_marker - marker turbo-performance value.
   * \param[in] val_span - span value.
   * \return The span-wise inflow tangential velocity.
   */
  inline su2double GetTangGridVelIn(unsigned short val_marker, unsigned short val_span) const {
    return TangGridVelIn[val_marker][val_span];
  }

  /*!
   * \brief Get the value of the outflow tangential velocity at each span.
   * \param[in] val_marker - marker turbo-performance value.
   * \param[in] val_span - span value.
   * \return The span-wise outflow tangential velocity.
   */
  inline su2double GetTangGridVelOut(unsigned short val_marker, unsigned short val_span) const {
    return TangGridVelOut[val_marker][val_span];
  }

  /*!
   * \brief Get the value of the inflow area at each span.
   * \param[in] val_marker - marker turbo-performance value.
   * \param[in] val_span - span value.
   * \return The span-wise inflow area.
   */
  inline su2double GetSpanAreaIn(unsigned short val_marker, unsigned short val_span) const {
    return SpanAreaIn[val_marker][val_span];
  }

  /*!
   * \brief Get the value of the outflow area at each span.
   * \param[in] val_marker - marker turbo-performance value.
   * \param[in] val_span - span value.
   * \return The span-wise outflow area.
   */
  inline su2double GetSpanAreaOut(unsigned short val_marker, unsigned short val_span) const {
    return SpanAreaOut[val_marker][val_span];
  }

  /*!
   * \brief Get the value of the inflow radius at each span.
   * \param[in] val_marker - marker turbo-performance value.
   * \param[in] val_span - span value.
   * \return The span-wise inflow radius.
   */
  inline su2double GetTurboRadiusIn(unsigned short val_marker, unsigned short val_span) const {
    return TurboRadiusIn[val_marker][val_span];
  }

  /*!
   * \brief Get the value of the outflow radius at each span.
   * \param[in] val_marker - marker turbo-performance value.
   * \param[in] val_span - span value.
   * \return The span-wise outflow radius.
   */
  inline su2double GetTurboRadiusOut(unsigned short val_marker, unsigned short val_span) const {
    return TurboRadiusOut[val_marker][val_span];
  }

  /*!
   * \brief Set the value of the inflow tangential velocity at each span.
   * \param[in] val_marker - marker turbo-performance value.
   * \param[in] val_span - span value.
   */
  inline void SetTangGridVelIn(su2double value, unsigned short val_marker, unsigned short val_span) {
    TangGridVelIn[val_marker][val_span] = value;
  }

  /*!
   * \brief Set the value of the outflow tangential velocity at each span.
   * \param[in] val_marker - marker turbo-performance value.
   * \param[in] val_span - span value.
   */
  inline void SetTangGridVelOut(su2double value, unsigned short val_marker, unsigned short val_span) {
    TangGridVelOut[val_marker][val_span] = value;
  }

  /*!
   * \brief Set the value of the inflow area at each span.
   * \param[in] val_marker - marker turbo-performance value.
   * \param[in] val_span - span value.
   */
  inline void SetSpanAreaIn(su2double value, unsigned short val_marker, unsigned short val_span) {
    SpanAreaIn[val_marker][val_span] = value;
  }

  /*!
   * \brief Set the value of the outflow area at each span.
   * \param[in] val_marker - marker turbo-performance value.
   * \param[in] val_span - span value.
   */
  inline void SetSpanAreaOut(su2double value, unsigned short val_marker, unsigned short val_span) {
    SpanAreaOut[val_marker][val_span] = value;
  }

  /*!
   * \brief Set the value of the inflow radius at each span.
   * \param[in] val_marker - marker turbo-performance value.
   * \param[in] val_span - span value.
   */
  inline void SetTurboRadiusIn(su2double value, unsigned short val_marker, unsigned short val_span) {
    TurboRadiusIn[val_marker][val_span] = value;
  }

  /*!
   * \brief Set the value of the outflow radius at each span.
   * \param[in] val_marker - marker turbo-performance value.
   * \param[in] val_span - span value.
   */
  inline void SetTurboRadiusOut(su2double value, unsigned short val_marker, unsigned short val_span) {
    TurboRadiusOut[val_marker][val_span] = value;
  }

  /*!
   * \brief A total number of vertex independently from the MPI partions.
   * \param[in] val_marker - marker value.
   * \param[in] val_span - span value.
   */
  inline unsigned long GetnTotVertexSpan(unsigned short val_marker, unsigned short val_span) const {
    return nTotVertexSpan[val_marker][val_span];
  }

  /*!
   * \brief min angular pitch independently from the MPI partions.
   * \param[in] val_marker - marker value.
   * \param[in] val_span - span value.
   */
  inline su2double GetMinAngularCoord(unsigned short val_marker, unsigned short val_span) const {
    return MinAngularCoord[val_marker][val_span];
  }

  /*!
   * \brief max angular pitch independently from the MPI partions.
   * \param[in] val_marker - marker value.
   * \param[in] val_span - span value.
   */
  inline su2double GetMaxAngularCoord(unsigned short val_marker, unsigned short val_span) const {
    return MaxAngularCoord[val_marker][val_span];
  }

  /*!
   * \brief min Relatice angular coord independently from the MPI partions.
   * \param[in] val_marker - marker value.
   * \param[in] val_span - span value.
   */
  inline su2double GetMinRelAngularCoord(unsigned short val_marker, unsigned short val_span) const {
    return MinRelAngularCoord[val_marker][val_span];
  }

  /*!
   * \brief Get the average grid velocity at a specific span for a given marker.
   * \param[in] val_marker - marker value.
   * \param[in] val_span - span value.
   */
  inline const su2double* GetAverageGridVel(unsigned short val_marker, unsigned short val_span) const {
    return AverageGridVel[val_marker][val_span];
  }

  /*!
   * \brief A virtual member.
   * \param config - Config
   */
  inline virtual void Check_Periodicity(CConfig* config) {}

  /*!
   * \brief Get the value of the customized temperature at a specified vertex on a specified marker.
   * \param[in] val_marker - Marker value
   * \param[in] val_vertex - Boundary vertex value
   */
  inline su2double GetCustomBoundaryTemperature(unsigned short val_marker, unsigned long val_vertex) const {
    return CustomBoundaryTemperature[val_marker][val_vertex];
  }

  /*!
   * \brief Set the value of the customized temperature at a specified vertex on a specified marker.
   * \param[in] val_marker - Marker value
   * \param[in] val_vertex - Boundary vertex value
   * \param[in] val_customBoundaryTemperature - Value of the temperature.
   */
  inline void SetCustomBoundaryTemperature(unsigned short val_marker, unsigned long val_vertex,
                                           su2double val_customBoundaryTemperature) {
    CustomBoundaryTemperature[val_marker][val_vertex] = val_customBoundaryTemperature;
  }

  /*!
   * \brief Get the value of the customized normal heat flux at a specified vertex on a specified marker.
   * \param[in] val_marker - Marker value
   * \param[in] val_vertex - Boundary vertex value
   */
  inline su2double GetCustomBoundaryHeatFlux(unsigned short val_marker, unsigned long val_vertex) const {
    return CustomBoundaryHeatFlux[val_marker][val_vertex];
  }

  /*!
   * \brief Set the value of the customized normal heat flux at a specified vertex on a specified marker.
   * \param[in] val_marker - Marker value
   * \param[in] val_vertex - Boundary vertex value
   * \param[in] val_customBoundaryHeatFlux - Value of the normal heat flux.
   */
  inline void SetCustomBoundaryHeatFlux(unsigned short val_marker, unsigned long val_vertex,
                                        su2double val_customBoundaryHeatFlux) {
    CustomBoundaryHeatFlux[val_marker][val_vertex] = val_customBoundaryHeatFlux;
  }

  /*!
   * \brief Filter values given at the element CG by performing a weighted average over a radial neighbourhood.
   * \param[in] filter_radius - Parameter defining the size of the neighbourhood.
   * \param[in] kernels - Kernel types and respective parameter, size of vector defines number of filter recursions.
   * \param[in] search_limit - Max degree of neighborhood considered for neighbor search, avoids excessive work in fine
   * regions. \param[in,out] values - On entry, the "raw" values, on exit, the filtered values.
   */
  void FilterValuesAtElementCG(const vector<su2double>& filter_radius,
                               const vector<pair<ENUM_FILTER_KERNEL, su2double>>& kernels,
                               const unsigned short search_limit, su2double* values) const;

  /*!
   * \brief Build the global (entire mesh!) adjacency matrix for the elements in compressed format.
   *        Used by FilterValuesAtElementCG to search for geometrically close neighbours.
   * \param[out] neighbour_start - i'th position stores the start position in "neighbour_idx" for the immediate
   *             neighbours of global element "i". Size nElemDomain+1
   * \param[out] neighbour_idx - Global index of the neighbours, mush be NULL on entry and free'd by calling function.
   */
  void GetGlobalElementAdjacencyMatrix(vector<unsigned long>& neighbour_start, long*& neighbour_idx) const;

  /*!
   * \brief Get the neighbours of the global element in the first position of "neighbours" that are within "radius" of
   * it. \param[in] iElem_global - Element of interest. \param[in] radius - Parameter defining the size of the
   * neighbourhood. \param[in] search_limit - Maximum "logical radius" to consider, limits cost in refined regions, use
   * 0 for unlimited. \param[in] neighbour_start - See GetGlobalElementAdjacencyMatrix. \param[in] neighbour_idx - See
   * GetGlobalElementAdjacencyMatrix. \param[in] cg_elem - Global element centroid coordinates in row major format
   * {x0,y0,x1,y1,...}. Size nDim*nElemDomain. \param[in,out] neighbours - The neighbours of iElem_global.
   * \param[in,out] is_neighbor - Working vector of size nElemGlobal, MUST be all false on entry (if so, on exit it will
   * be the same). \return true if the search was successful, i.e. not limited.
   */
  bool GetRadialNeighbourhood(const unsigned long iElem_global, const passivedouble radius, size_t search_limit,
                              const vector<unsigned long>& neighbour_start, const long* neighbour_idx,
                              const su2double* cg_elem, vector<long>& neighbours, vector<bool>& is_neighbor) const;

  /*!
   * \brief Compute and store the volume of the primal elements.
   */
  void SetElemVolume();

  /*!
   * \brief Set the multigrid index for the current geometry object.
   * \param[in] val_iMesh - Multigrid index for current geometry object.
   */
  inline void SetMGLevel(unsigned short val_iMesh) { MGLevel = val_iMesh; }

  /*!
   * \brief Get the multigrid index for the current geometry object.
   * \return Multigrid index for current geometry object.
   */
  inline unsigned short GetMGLevel() const { return MGLevel; }

  /*!
   * \brief A virtual member.
   * \param config - Config
   */
  inline virtual void ComputeMeshQualityStatistics(const CConfig* config) {}

  /*!
   * \brief Color multigrid levels for visualization.
   * \param nMGLevels - Number of levels
   * \param geometry - The levels
   */
  void ColorMGLevels(unsigned short nMGLevels, const CGeometry* const* geometry);

  /*!
   * \brief Get the sparse pattern of "type" with given level of fill.
   * \note This method builds the pattern if that has not been done yet.
   * \param[in] type - Finite volume or finite element.
   * \param[in] fillLvl - Level of fill of the pattern.
   * \return Reference to the sparse pattern.
   */
  const CCompressedSparsePatternUL& GetSparsePattern(ConnectivityType type, unsigned long fillLvl = 0);

  /*!
   * \brief Get the edge to sparse pattern map.
   * \note This method builds the map and required pattern (0-fill FVM) if that has not been done yet.
   * \return Reference to the map.
   */
  const CEdgeToNonZeroMapUL& GetEdgeToSparsePatternMap();

  /*!
   * \brief Get the transpose of the (main, i.e 0 fill) sparse pattern (e.g. CSR becomes CSC).
   * \param[in] type - Finite volume or finite element.
   * \return Reference to the map.
   */
  const su2vector<unsigned long>& GetTransposeSparsePatternMap(ConnectivityType type);

  /*!
   * \brief Get the edge coloring.
   * \note This method computes the coloring if that has not been done yet.
   * \param[out] efficiency - optional output of the coloring efficiency.
   * \return Reference to the coloring.
   */
  const CCompressedSparsePatternUL& GetEdgeColoring(su2double* efficiency = nullptr);

  /*!
   * \brief Force the natural (sequential) edge coloring.
   */
  void SetNaturalEdgeColoring();

  /*!
   * \brief Get the group size used in edge coloring.
   * \return Group size.
   */
  inline unsigned long GetEdgeColorGroupSize() const { return edgeColorGroupSize; }

  /*!
   * \brief Get the element coloring.
   * \note This method computes the coloring if that has not been done yet.
   * \param[out] efficiency - optional output of the coloring efficiency.
   * \return Reference to the coloring.
   */
  const CCompressedSparsePatternUL& GetElementColoring(su2double* efficiency = nullptr);

  /*!
   * \brief Force the natural (sequential) element coloring.
   */
  void SetNaturalElementColoring();

  /*!
   * \brief Get the group size used in element coloring.
   * \return Group size.
   */
  inline unsigned long GetElementColorGroupSize() const { return elemColorGroupSize; }

  /*!
   * \brief Get the linelet definition, this function computes the linelets if that has not been done yet.
   */
  const CLineletInfo& GetLineletInfo(const CConfig* config) const;

  /*!
   * \brief Compute an ADT including the coordinates of all viscous markers
   * \param[in] config - Definition of the particular problem.
   * \return pointer to the ADT
   */
  virtual std::unique_ptr<CADTElemClass> ComputeViscousWallADT(const CConfig* config) const { return nullptr; }

  /*!
   * \brief Reduce the wall distance based on an previously constructed ADT.
   * \details The ADT might belong to another zone, giving rise to lower wall distances
   * than those already stored.
   * \param[in] WallADT - The ADT to reduce the wall distance
   * \param[in] config - Config of this geometry (not the ADT zone's geometry)
   * \param[in] iZone - Zone whose markers made the ADT
   */
  virtual void SetWallDistance(CADTElemClass* WallADT, const CConfig* config,
                               unsigned short iZone = numeric_limits<unsigned short>::max()) {}

  /*!
   * \brief Set wall distances a specific value
   *  \param[in] val - new value for the wall distance at all points.
   */
  virtual void SetWallDistance(su2double val) {}

  /*!
   * \brief Compute the distances to the closest vertex on viscous walls over the entire domain
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] geometry_container - Geometrical definition of the problem.
   */
  static void ComputeWallDistance(const CConfig* const* config_container, CGeometry**** geometry_container);

  /*!
   * \brief Set the amount of nonconvex elements in the mesh.
   * \param[in] nonconvex_elems - amount of nonconvex elements in the mesh
   */
  void SetnNonconvexElements(unsigned long nonconvex_elems) { nNonconvexElements = nonconvex_elems; }

  /*!
   * \brief Get the amount of nonconvex elements in the mesh.
   * \param[out] nNonconvexElements- amount of nonconvex elements in the mesh
   */
  unsigned long GetnNonconvexElements() const { return nNonconvexElements; }

  /*!
   * \brief For streamwise periodicity, find & store a unique reference node on the designated periodic inlet.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void FindUniqueNode_PeriodicBound(const CConfig* config) {}

  /*!
   * \brief Get a pointer to the reference node coordinate vector.
   * \return A pointer to the reference node coordinate vector.
   */
  inline virtual const su2double* GetStreamwise_Periodic_RefNode() const { return nullptr; }
};
