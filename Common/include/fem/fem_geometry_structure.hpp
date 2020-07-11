/*!
 * \file fem_geometry_structure.hpp
 * \brief Headers of the main subroutines for creating the geometrical structure for the FEM solver.
 *        The subroutines and functions are in the <i>fem_geometry_structure.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../geometry/CGeometry.hpp"
#include "fem_standard_element.hpp"
#include "../wall_model.hpp"
#include "../blas_structure.hpp"

using namespace std;

/*!
 * \class CLong3T
 * \brief Help class used to store three longs as one entity.
 * \version 7.0.6 "Blackbird"
 */
struct CLong3T {
  long long0 = 0;  /*!< \brief First long to store in this class. */
  long long1 = 0;  /*!< \brief Second long to store in this class. */
  long long2 = 0;  /*!< \brief Third long to store in this class. */

  CLong3T() = default;

  CLong3T(const long a, const long b, const long c) {long0 = a; long1 = b; long2 = c;}

  bool operator<(const CLong3T &other) const;
};

/*!
 * \class CVolumeElementFEM
 * \brief Class to store a volume element for the FEM solver.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
class CVolumeElementFEM {
public:

  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
  unsigned short nPolyGrid;    /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short nPolySol;     /*!< \brief Polynomial degree for the solution of the element. */
  unsigned short nDOFsGrid;    /*!< \brief Number of DOFs for the geometry of the element. */
  unsigned short nDOFsSol;     /*!< \brief Number of DOFs for the solution of the element. */

  unsigned short indStandardElement; /*!< \brief Index in the vector of standard elements. */

  unsigned long offsetDOFsSolGlobal; /*!< \brief Global offset of the solution DOFs of this element. */
  unsigned long offsetDOFsSolLocal;  /*!< \brief Local offset of the solution DOFs of this element. */

  vector<su2double> coorSolDOFs;            /*!< \brief The coordinates of the solution DOFs of this element. */

  /*!
   * \brief Get all the corner points of all the faces of this element. It must be made sure
            that the numbering of the faces is identical to the numbering used for the
            standard elements.
   * \param[out] nFaces         - Number of faces of this element.
   * \param[out] nPointsPerFace - Number of corner points for each of the faces.
   * \param[out] faceConn       - Global IDs of the corner points of the faces.
   */
  void GetCornerPointsAllFaces(unsigned short &numFaces,
                               unsigned short nPointsPerFace[],
                               unsigned long  faceConn[6][4]);
};

/*!
 * \class CPointFEM
 * \brief Class to a point for the FEM solver.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
struct CPointFEM {
  unsigned long globalID;    /*!< \brief The global ID of this point in the grid. */
  short periodIndexToDonor;  /*!< \brief The index of the periodic transformation to the donor
                                  element. Only for halo elements. A -1 indicates no periodic
                                  transformation. */
  su2double coor[3] = {0.0}; /*!< \brief Array with the coordinates of the node. */

  /*!
   * \brief Less than operator of the class. Needed for the sorting.
   */
  bool operator<(const CPointFEM &other) const;

  /*!
   * \brief Equal operator of the class. Needed for the removal of double entities.
   */
  bool operator==(const CPointFEM &other) const;

};

/*!
 * \class CInternalFaceElementFEM
 * \brief Class to store an internal face for the FEM solver.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
struct CInternalFaceElementFEM {

};

/*!
 * \class CSurfaceElementFEM
 * \brief Class to store a surface element for the FEM solver.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
struct CSurfaceElementFEM {

  unsigned short indStandardElement; /*!< \brief Index in the vector of standard elements. */

  vector<unsigned long> DOFsSolFace;    /*!< \brief Vector of the solution DOFs of the face. */

  /*!
   *  \brief Function, which determines the corner points of this surface element.
   *  \param[out] nPointsPerFace - Number of corner points of the face.
   *  \param[out] faceConn       - The corner points of the face.
   */
  void GetCornerPointsFace(unsigned short &nPointsPerFace,
                           unsigned long  faceConn[]);
};

/*!
 * \class CBoundaryFEM
 * \brief Class to store a boundary for the FEM solver.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
struct CBoundaryFEM {
  string markerTag;  /*!< \brief Marker tag of this boundary. */

  bool periodicBoundary = false;     /*!< \brief Whether or not this boundary is a periodic boundary. */
  bool haloInfoNeededForBC = false;  /*!< \brief Whether or not information of halo elements
                                                 is needed to impose the boundary conditions. */

  vector<unsigned long> nSurfElem; /*!< \brief Number of surface elements per time level,
                                               cumulative storage format. */

  vector<CSurfaceElementFEM> surfElem; /*!< \brief Vector of the local surface elements. */

  CWallModel *wallModel = nullptr;     /*!< \brief Wall model for LES. */

  ~CBoundaryFEM(void) { delete wallModel; }
};

/*!
 * \class CMeshFEM
 * \brief Base class for the FEM solver.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
class CMeshFEM: public CGeometry {
protected:
  unsigned long nVolElemTot{0};    /*!< \brief Total number of local volume elements, including halos. */
  unsigned long nVolElemOwned{0};  /*!< \brief Number of owned local volume elements. */

  vector<CVolumeElementFEM> volElem; /*!< \brief Vector of the local volume elements, including halos. */

  vector<CBoundaryFEM> boundaries;   /*!< \brief Vector of the boundaries of the FEM mesh. */

  vector<int> ranksRecv;             /*!< \brief Vector of ranks, from which this rank will receive halo
                                                 information. Self communication is included. */
  vector<int> ranksSend;             /*!< \brief Vector of ranks, to which this rank will send halo
                                                 information. Self communication is included. */

  vector<vector<unsigned long> > entitiesSend;    /*!< \brief Vector of vector, which contains the entities that
                                                              must be sent. Self communication is included. For DG
                                                              an entitity is an element, for regular FEM an entity
                                                              is a DOF. */
  vector<vector<unsigned long> > entitiesRecv; /*!< \brief Vector of vector, which contains the entities that
                                                           must be received. Self communication is included. For DG
                                                           an entity is an element, for regular FEM an entity
                                                           is a DOF. */

  vector<CFEMStandardBoundaryFace> standardBoundaryFacesSol;  /*!< \brief Vector that contains the standard boundary
                                                                          faces used for the solution of the DG solver. */
  vector<CFEMStandardBoundaryFace> standardBoundaryFacesGrid; /*!< \brief Vector that contains the standard boundary
                                                                          faces used for the geometry of the DG solver. */

  CBlasStructure *blasFunctions{nullptr}; /*!< \brief  Pointer to the object to carry out the BLAS functionalities. */

public:
  /*!
  * \brief Constructor of the class.
  */
 CMeshFEM(void) : CGeometry() { }

  /*!
  * \overload
  * \brief Redistributes the grid over the ranks and creates the halo layer.
  * \param[in] geometry - The linear distributed grid that must be redistributed.
  * \param[in] config   - Definition of the particular problem.
  */
  CMeshFEM(CGeometry *geometry, CConfig *config);

  /*!
  * \brief Destructor of the class.
  */
  ~CMeshFEM(void) override { delete blasFunctions; }

  /*!
  * \brief Function, which makes available the boundaries of the local FEM mesh.
  * \return  Pointer to the boundaries of the local FEM mesh.
  */
  inline CBoundaryFEM* GetBoundaries(void) {return boundaries.data();}

  /*!
  * \brief Function, which makes available the number of owned volume elements in the local FEM mesh.
  * \return  Number of owned volume elements of the local FEM mesh.
  */
  inline unsigned long GetNVolElemOwned(void) const {return nVolElemOwned;}

  /*!
  * \brief Function, which makes available the total number of volume elements in the local FEM mesh.
  * \return  Total number of volume elements of the local FEM mesh.
  */
  inline unsigned long GetNVolElemTot(void) const {return nVolElemTot;}

  /*!
  * \brief Function, which makes available the volume elements in the local FEM mesh.
  * \return  Pointer to the volume elements of the local FEM mesh.
  */
  inline CVolumeElementFEM* GetVolElem(void) {return volElem.data();}

  /*!
  * \brief Function, which makes available the number of standard boundary faces of the solution.
  * \return  Number of standard boundary faces of the solution.
  */
  inline unsigned short GetNStandardBoundaryFacesSol(void) {return standardBoundaryFacesSol.size();}

  /*!
  * \brief Function, which makes available the standard boundary faces of the solution.
  * \return  Pointer to the standard boundary faces of the solution.
  */
  inline CFEMStandardBoundaryFace* GetStandardBoundaryFacesSol(void) {return standardBoundaryFacesSol.data();}

  /*!
  * \brief Function, which makes available the vector of receive ranks as
           a const reference.
  * \return  Const reference to the vector of ranks.
  */
  inline const vector<int>& GetRanksRecv(void) const {return ranksRecv;}

  /*!
  * \brief Function, which makes available the vector of send ranks as
           a const reference.
  * \return  Const reference to the vector of ranks.
  */
  inline const vector<int>& GetRanksSend(void) const {return ranksSend;}

  /*!
  * \brief Function, which makes available the vector of vectors containing the receive
           entities as a const reference.
  * \return  Const reference to the vector of vectors of receive entities.
  */
  inline const vector<vector<unsigned long> >& GetEntitiesRecv(void) const {return entitiesRecv;}

  /*!
  * \brief Function, which makes available the vector of vectors containing the send
           entities as a const reference.
  * \return  Const reference to the vector of vectors of send entities.
  */
  inline const vector<vector<unsigned long> >& GetEntitiesSend(void) const {return entitiesSend;}

  /*!
  * \brief Compute surface area (positive z-direction) for force coefficient non-dimensionalization.
  * \param[in] config - Definition of the particular problem.
  */
  void SetPositive_ZArea(CConfig *config) override;

protected:
  /*!
  * \brief Function, which computes the gradients of the parametric coordinates
           w.r.t. the Cartesian coordinates in the integration points of a face,
           i.e. drdx, drdy, dsdx, etc.
  * \param[in]  nIntegration   - Number of integration points on the face.
  * \param[in]  nDOFs          - Number of DOFs of the grid associated with the
                                 neighboring element.
  * \param[in]  matDerBasisInt - Matrix, which contains the derivatives of the
                                 basis functions w.r.t. the parametric
                                 coordinates r, s and t in the integration points.
  * \param[in]  DOFs           - The DOFs of the grid associated with the element.
  * \param[out] derivCoor      - Storage for the derivatives of the coordinates.
  *  \param[in] config         - Definition of the particular problem.
  */
  void ComputeGradientsCoordinatesFace(const unsigned short nIntegration,
                                       const unsigned short nDOFs,
                                       const su2double      *matDerBasisInt,
                                       const unsigned long  *DOFs,
                                       su2double            *derivCoor,
                                       CConfig              *config);
  /*!
  * \brief Function, which computes the gradients of the Cartesian coordinates
           w.r.t. the parametric coordinates in the given set of integration
           points, i.e. dxdr, dydr, etc.
  * \param[in]  nIntegration   - Number of integration points.
  * \param[in]  nDOFs          - Number of DOFs of the grid associated with the
                                 element.
  * \param[in]  matDerBasisInt - Matrix, which contains the derivatives of the
                                 basis functions w.r.t. the parametric
                                 coordinates r, s and t in the integration points.
  * \param[in]  DOFs           - The DOFs of the grid associated with the element.
  * \param[out] derivCoor    - Storage for the derivatives of the coordinates.
  * \param[in] config        - Definition of the particular problem.
  */
  void ComputeGradientsCoorWRTParam(const unsigned short nIntegration,
                                    const unsigned short nDOFs,
                                    const su2double      *matDerBasisInt,
                                    const unsigned long  *DOFs,
                                    su2double            *derivCoor,
                                    CConfig              *config);
  /*!
  * \brief Function, which computes the information of the normals in the
           integration points of a face.
  * \param[in]  nIntegration - Number of integration points on the face.
  * \param[in]  nDOFs        - Number of DOFs of the grid associated with the face.
  * \param[in]  dr           - r-derivatives of the basis functions of the face.
  * \param[in]  ds           - s-derivatives of the basis functions of the face.
                               Only for 3D computations.
  * \param[in]  DOFs         - The DOFs of the grid associated with the face.
  * \param[out] normals      - Storage for the normal information to be computed.
  */
  void ComputeNormalsFace(const unsigned short nIntegration,
                          const unsigned short nDOFs,
                          const su2double      *dr,
                          const su2double      *ds,
                          const unsigned long  *DOFs,
                          su2double            *normals);

  /*!
  * \brief Function, which computes the metric terms of the faces of a
           physical boundary.
  * \param[inout] boundary - Boundary for whose faces the boundary metric
                             terms must be computed.
  * \param[in]    config   - Definition of the particular problem.
  */
  void MetricTermsBoundaryFaces(CBoundaryFEM *boundary,
                                CConfig      *config);
};

/*!
 * \class CMeshFEM_DG
 * \brief Class which contains all the variables for the DG FEM solver.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
class CMeshFEM_DG: public CMeshFEM {

public:
  /*!
   * \brief Constructor of the class.
   */
  CMeshFEM_DG(void) : CMeshFEM() {}

  /*!
  * \overload
  * \brief Redistributes the grid over the ranks and creates the halo layer.
  * \param[in] geometry - The linear distributed grid that must be redistributed.
  * \param[in] config   - Definition of the particular problem.
  */
  CMeshFEM_DG(CGeometry *geometry, CConfig *config);

 /*!
  * \brief Function to compute the coordinates of the integration points.
  */
  void CoordinatesIntegrationPoints(void);

 /*!
  * \brief Function to compute the coordinates of solution DOFs.
  */
  void CoordinatesSolDOFs(void);

 /*!
  * \brief Function to create the faces used in the DG formulation.
  * \param[in] config - Definition of the particular problem.
  */
  void CreateFaces(CConfig *config);

 /*!
  * \brief Function to create the standard volume elements.
  * \param[in] config - Definition of the particular problem.
  */
  void CreateStandardVolumeElements(CConfig *config);

 /*!
  * \brief Function to compute the grid velocities for static problems.
  * \param[in] config             - Definition of the particular problem.
  * \param[in] Kind_Grid_Movement - The type of prescribed grid motion.
  * \param[in] iZone              - The currently active zone number.
  */
  void InitStaticMeshMovement(CConfig              *config,
                              const unsigned short Kind_Grid_Movement,
                              const unsigned short iZone);

 /*!
  * \brief Function, which makes available the standard volume elements of the solution.
  * \return  Pointer to the standard volume elements of the solution.
  */
  inline CFEMStandardElement* GetStandardElementsSol(void) {return NULL;}

 /*!
  * \brief Function, which computes a length scale of the volume elements.
  * \param[in] config - Definition of the particular problem.
  */
  void LengthScaleVolumeElements(void);

 /*!
  * \brief Function, which computes the metric terms of the surface
           elements, both internal faces and physical boundary faces.
  * \param[in] config - Definition of the particular problem.
  */
  void MetricTermsSurfaceElements(CConfig *config);

 /*!
  * \brief Function, which computes the metric terms of the
           volume elements.
  * \param[in] config - Definition of the particular problem.
  */
  void MetricTermsVolumeElements(CConfig *config);

 /*!
  * \brief Set the send receive boundaries of the grid.
  * \param[in] config - Definition of the particular problem.
  */
  void SetSendReceive(CConfig *config) override;

  /*!
   * \brief Set the local index that correspond with the global numbering index.
   */
  void SetGlobal_to_Local_Point() override;

  /*!
   * \brief Function, which carries out the preprocessing tasks when wall functions are used.
   * \param[in] config - Definition of the particular problem.
   */
  void WallFunctionPreprocessing(CConfig *config);

protected:
 /*!
  * \brief Function, which computes the correct sequence of the connectivities
           of a face, such that it matches the sequence of the given corner points.
  * \param[in]  VTK_TypeFace       - Type of the face using the VTK convention.
  * \param[in]  cornerPointsFace   - The corner points of the face in the desired
                                     sequence.
  * \param[in]  VTK_TypeElem       - Type of the element using the VTK convention.
  * \param[in]  nPolyGrid          - Polynomial degree used in the grid definition
                                     for the face and the element.
  * \param[in]  elemNodeIDsGrid    - The node IDs of the grid DOFs of the element,
                                     i.e. the element connectivity.
  * \param[in]  nPolyConn          - Polynomial degree of the connectivities to
                                     be modified.
  * \param[in]  connElem           - Connectivity of the adjacent volume element.
  * \param[out] swapFaceInElement  - Whether or not the connectivity of the face must
                                     be swapped compared to the face of the corresponding
                                     standard element. Only relevant for triangular faces
                                     of a pyramid and quadrilateral faces of a prism.
                                     corresponds to the top point of the adjacent pyramid.
  * \param[out] modConnFace        - Connectivity of the face after the renumbering.
  * \param[out] modConnElem        - Connectivity of the element after the renumbering.
                                     This renumbering is such that the face corresponds
                                     to the appropriate face of the element used in the
                                     standard faces and also the corner points match.
  */
  void CreateConnectivitiesFace(const unsigned short        VTK_TypeFace,
                                const unsigned long         *cornerPointsFace,
                                const unsigned short        VTK_TypeElem,
                                const unsigned short        nPolyGrid,
                                const vector<unsigned long> &elemNodeIDsGrid,
                                const unsigned short        nPolyConn,
                                const unsigned long         *connElem,
                                bool                        &swapFaceInElement,
                                unsigned long               *modConnFace,
                                unsigned long               *modConnElem);

  /*!
  * \brief Function, which computes the correct sequence of the connectivities
           of a line adjacent to a quadrilateral, such that the line is face 0
           of the quadrilateral and it matches the sequence of the given corner points.
  * \param[in]  cornerPointsLine - The corner points of the line in the desired
                                   sequence.
  * \param[in]  nPolyGrid        - Polynomial degree used in the grid definition
                                   for the line and the quadrilateral.
  * \param[in]  quadNodeIDsGrid  - The node IDs of the grid DOFs of the quadrilateral,
                                   i.e. the element connectivity.
  * \param[in]  nPolyConn        - Polynomial degree of the connectivities to
                                   be modified.
  * \param[in]  connQuad         - Connectivity of the adjacent quadrilateral.
  * \param[out] modConnLine      - Connectivity of the line after the renumbering.
  * \param[out] modConnQuad      - Connectivity of the quadrilateral after the
                                   renumbering. This renumbering is such that the
                                   line corresponds to face 0 of the quadrilateral.
  */
  void CreateConnectivitiesLineAdjacentQuadrilateral(
                                  const unsigned long         *cornerPointsLine,
                                  const unsigned short        nPolyGrid,
                                  const vector<unsigned long> &quadNodeIDsGrid,
                                  const unsigned short        nPolyConn,
                                  const unsigned long         *connQuad,
                                  unsigned long               *modConnLine,
                                  unsigned long               *modConnQuad);

  /*!
  * \brief Function, which computes the correct sequence of the connectivities
           of a line adjacent to a triangle, such that the line is face 0
           of the triangle and it matches the sequence of the given corner points.
  * \param[in]  cornerPointsLine - The corner points of the line in the desired
                                   sequence.
  * \param[in]  nPolyGrid        - Polynomial degree used in the grid definition
                                   for the line and the triangle.
  * \param[in]  triaNodeIDsGrid  - The node IDs of the grid DOFs of the triangle,
                                   i.e. the element connectivity.
  * \param[in]  nPolyConn        - Polynomial degree of the connectivities to
                                   be modified.
  * \param[in]  connTria         - Connectivity of the adjacent triangle.
  * \param[out] modConnLine      - Connectivity of the line after the renumbering.
  * \param[out] modConnTria      - Connectivity of the triangle after the
                                   renumbering. This renumbering is such that the
                                   line corresponds to face 0 of the triangle.
  */
  void CreateConnectivitiesLineAdjacentTriangle(
                                  const unsigned long         *cornerPointsLine,
                                  const unsigned short        nPolyGrid,
                                  const vector<unsigned long> &triaNodeIDsGrid,
                                  const unsigned short        nPolyConn,
                                  const unsigned long         *connTria,
                                  unsigned long               *modConnLine,
                                  unsigned long               *modConnTria);

  /*!
  * \brief Function, which computes the correct sequence of the connectivities
           of a quad adjacent to a hexahedron, such that the quad is face 0
           of the hexahedron and it matches the sequence of the given corner points.
  * \param[in]  cornerPointsQuad - The corner points of the quad in the desired
                                   sequence.
  * \param[in]  nPolyGrid        - Polynomial degree used in the grid definition
                                   for the quadrilateral and the hexahedron.
  * \param[in]  hexaNodeIDsGrid  - The node IDs of the grid DOFs of the
                                   hexahedron, i.e. the element connectivity.
  * \param[in]  nPolyConn        - Polynomial degree of the connectivities to
                                   be modified.
  * \param[in]  connHexa         - Connectivity of the adjacent hexahedron.
  * \param[out] modConnQuad      - Connectivity of the quad after the renumbering.
  * \param[out] modConnHexa      - Connectivity of the hexahedron after the
                                   renumbering. This renumbering is such that the
                                   quad corresponds to face 0 of the hexahedron.
  */
  void CreateConnectivitiesQuadrilateralAdjacentHexahedron(
                                  const unsigned long         *cornerPointsQuad,
                                  const unsigned short        nPolyGrid,
                                  const vector<unsigned long> &hexaNodeIDsGrid,
                                  const unsigned short        nPolyConn,
                                  const unsigned long         *connHexa,
                                  unsigned long               *modConnQuad,
                                  unsigned long               *modConnHexa);

  /*!
  * \brief Function, which computes the correct sequence of the connectivities
           of a quad adjacent to a prism, such that the quad is face 2
           of the prism and it matches the sequence of the given corner points.
  * \param[in]  cornerPointsQuad  - The corner points of the quad in the desired
                                    sequence.
  * \param[in]  nPolyGrid         - Polynomial degree used in the grid definition
                                    for the quadrilateral and the prism.
  * \param[in]  prismNodeIDsGrid  - The node IDs of the grid DOFs of the prism,
                                    i.e. the element connectivity.
  * \param[in]  nPolyConn         - Polynomial degree of the connectivities to
                                    be modified.
  * \param[in]  connPrism         - Connectivity of the adjacent prism.
  * \param[out] swapFaceInElement - Whether or not the connectivity of the face must
                                    be swapped compared to the face of the corresponding
                                    standard element.
  * \param[out] modConnQuad       - Connectivity of the quad after the renumbering.
  * \param[out] modConnPrism      - Connectivity of the prism after the
                                    renumbering. This renumbering is such that the
                                    quad corresponds to face 3 of the prism.
  */
  void CreateConnectivitiesQuadrilateralAdjacentPrism(
                                  const unsigned long         *cornerPointsQuad,
                                  const unsigned short        nPolyGrid,
                                  const vector<unsigned long> &prismNodeIDsGrid,
                                  const unsigned short        nPolyConn,
                                  const unsigned long         *connPrism,
                                  bool                        &swapFaceInElement,
                                  unsigned long               *modConnQuad,
                                  unsigned long               *modConnPrism);

  /*!
  * \brief Function, which computes the correct sequence of the connectivities
           of a quad adjacent to a pyramid, such that the quad is face 0
           of the pyramid and it matches the sequence of the given corner points.
  * \param[in]  cornerPointsQuad - The corner points of the quad in the desired
                                   sequence.
  * \param[in]  nPolyGrid        - Polynomial degree used in the grid definition
                                   for the quadrilateral and the pyramid.
  * \param[in]  pyraNodeIDsGrid  - The node IDs of the grid DOFs of the pyramid,
                                   i.e. the element connectivity.
  * \param[in]  nPolyConn        - Polynomial degree of the connectivities to
                                   be modified.
  * \param[in]  connPyra         - Connectivity of the adjacent pyramid.
  * \param[out] modConnQuad      - Connectivity of the quad after the renumbering.
  * \param[out] modConnPyra      - Connectivity of the pyramid after the
                                   renumbering. This renumbering is such that the
                                   quad corresponds to face 0 of the pyramid.
  */
  void CreateConnectivitiesQuadrilateralAdjacentPyramid(
                                  const unsigned long         *cornerPointsQuad,
                                  const unsigned short        nPolyGrid,
                                  const vector<unsigned long> &pyraNodeIDsGrid,
                                  const unsigned short        nPolyConn,
                                  const unsigned long         *connPyra,
                                  unsigned long               *modConnQuad,
                                  unsigned long               *modConnPyra);

  /*!
  * \brief Function, which computes the correct sequence of the connectivities
           of a triangle adjacent to a prism, such that the triangle is face 0
           of the prism and it matches the sequence of the given corner points.
  * \param[in]  cornerPointsTria - The corner points of the triangle in the
                                   desired sequence.
  * \param[in]  nPolyGrid        - Polynomial degree used in the grid definition
                                   for the triangle and the prism.
  * \param[in]  prismNodeIDsGrid - The node IDs of the grid DOFs of the prism,
                                   i.e. the element connectivity.
  * \param[in]  nPolyConn        - Polynomial degree of the connectivities to
                                   be modified.
  * \param[in]  connPrism        - Connectivity of the adjacent prism.
  * \param[out] modConnTria      - Connectivity of the triangle after the renumbering.
  * \param[out] modConnPrism     - Connectivity of the prism after the
                                   renumbering. This renumbering is such that the
                                   triangle corresponds to face 0 of the prism.
  */
  void CreateConnectivitiesTriangleAdjacentPrism(
                                  const unsigned long         *cornerPointsTria,
                                  const unsigned short        nPolyGrid,
                                  const vector<unsigned long> &prismNodeIDsGrid,
                                  const unsigned short        nPolyConn,
                                  const unsigned long         *connPrism,
                                  unsigned long               *modConnTria,
                                  unsigned long               *modConnPrism);

  /*!
  * \brief Function, which computes the correct sequence of the connectivities
           of a triangle adjacent to a pyramid, such that the triangle is face 1
           of the pyramid and it matches the sequence of the given corner points.
  * \param[in]  cornerPointsTria   - The corner points of the triangle in the
                                     desired sequence.
  * \param[in]  nPolyGrid          - Polynomial degree used in the grid definition
                                     for the triangle and the pyramid.
  * \param[in]  pyraNodeIDsGrid    - The node IDs of the grid DOFs of the pyramid,
                                     i.e. the element connectivity.
  * \param[in]  nPolyConn          - Polynomial degree of the connectivities to
                                     be modified.
  * \param[in]  connPyra           - Connectivity of the adjacent pyramid.
  * \param[out] swapFaceInElement  - Whether or not the connectivity of the face must
                                     be swapped compared to the face of the corresponding
                                     standard element.
  * \param[out] modConnTria        - Connectivity of the triangle after the renumbering.
  * \param[out] modConnPyra        - Connectivity of the pyramid after the
                                     renumbering. This renumbering is such that the
                                     triangle corresponds to face 3 of the pyramid.
  */
  void CreateConnectivitiesTriangleAdjacentPyramid(
                                  const unsigned long         *cornerPointsTria,
                                  const unsigned short        nPolyGrid,
                                  const vector<unsigned long> &pyraNodeIDsGrid,
                                  const unsigned short        nPolyConn,
                                  const unsigned long         *connPyra,
                                  bool                        &swapFaceInElement,
                                  unsigned long               *modConnTria,
                                  unsigned long               *modConnPyra);

  /*!
  * \brief Function, which computes the correct sequence of the connectivities
           of a triangle adjacent to a tetrahedron, such that the triangle is face 0
           of the tetrahedron and it matches the sequence of the given corner points.
  * \param[in]  cornerPointsTria - The corner points of the triangle in the
                                   desired sequence.
  * \param[in]  nPolyGrid        - Polynomial degree used in the grid definition
                                   for the triangle and the tetrahedron.
  * \param[in]  tetNodeIDsGrid   - The node IDs of the grid DOFs of the
                                   tetrahedron, i.e. the element connectivity.
  * \param[in]  nPolyConn        - Polynomial degree of the connectivities to
                                   be modified.
  * \param[in]  connTet          - Connectivity of the adjacent tetrahedron.
  * \param[out] modConnTria      - Connectivity of the triangle after the renumbering.
  * \param[out] modConnTet       - Connectivity of the tetrahedron after the
                                   renumbering. This renumbering is such that the
                                   triangle corresponds to face 0 of the tetrahedron.
  */
  void CreateConnectivitiesTriangleAdjacentTetrahedron(
                                  const unsigned long         *cornerPointsTria,
                                  const unsigned short        nPolyGrid,
                                  const vector<unsigned long> &tetNodeIDsGrid,
                                  const unsigned short        nPolyConn,
                                  const unsigned long         *connTet,
                                  unsigned long               *modConnTria,
                                  unsigned long               *modConnTet);

 /*!
  * \brief Function, which computes the parametric coordinates of the given
           Cartesian coordinates inside the given parent element.
  * \param[in]  coor           - Cartesian coordinates for which the parametric
                                 coordinates must be determined.
  * \param[in]  parElem        - The high order parent element which contains
                                 the point.
  * \param[in]  subElem        - Low order sub element inside the parent element
                                 which contains the point.
  * \param[in]  weightsSubElem - Interpolation weights inside subElem for the
                                 coordinates. Used for an initial guess.
  * \param[out] parCoor        - Parametric coordinates inside the high order
                                 parent element for the given coordinates.
                                 These parametric coordinates must be computed.
  */
  void HighOrderContainmentSearch(const su2double      *coor,
                                  const unsigned long  parElem,
                                  const unsigned short subElem,
                                  const su2double      *weightsSubElem,
                                  su2double            *parCoor);

  /*!
  * \brief Function, which computes the metric terms for internal
           matching faces.
  * \param[in] config - Definition of the particular problem.
  */
  void MetricTermsMatchingFaces(CConfig *config);

  /*!
  * \brief Function, which computes the time coefficients for the ADER-DG predictor step.
  * \param[in] config - Definition of the particular problem.
  */
  void TimeCoefficientsPredictorADER_DG(CConfig *config);

  /*!
  * \brief Function, which computes the volume metric terms for the given
           entities from the gradients of the coordinates.
  * \param[in]  nEntities    - Number of entities for which the metric terms
                               must be computed.
  * \param[in]  gradCoor     - The gradients of the coordinates (w.r.t. the
                               parametric coordinates) from which the metric
                               terms must be computed.
  * \param[out] metricTerms  - Vector in which the metric terms must be stored.
  */
  void VolumeMetricTermsFromCoorGradients(const unsigned short nEntities,
                                          const su2double      *gradCoor,
                                          vector<su2double>    &metricTerms);

  /*!
   * \brief Compute an ADT including the coordinates of all viscous markers
   * \param[in] config - Definition of the particular problem.
   * \return pointer to the ADT
   */
  std::unique_ptr<CADTElemClass> ComputeViscousWallADT(const CConfig *config) const override;

  /*!
   * \brief Set wall distances a specific value
   */
  void SetWallDistance(su2double val) override;

  /*!
   * \brief Set the wall distance based on an previously constructed ADT
   * \param[in] config - Definition of the particular problem.
   * \param[in] WallADT - The ADT to compute the wall distance
   */
  void SetWallDistance(const CConfig *config, CADTElemClass* WallADT) override;
};

/*!
 * \class CDummyMeshFEM_DG
 * \brief Class for defining a DG geometry that does not contain any points/elements.
 *        Can be used for initializing other classes that depend on the geometry without
 *        going through the time-consuming mesh initialization and paritioning.
 * \author T. Albring
 */
class CDummyMeshFEM_DG : public CMeshFEM_DG {

public:
  /*!
   * \brief Constructor of the class
   * \param[in] config - Definition of the particular problem.
   */
  CDummyMeshFEM_DG(CConfig *config);

};
