/*!
 * \file fem_geometry_structure.hpp
 * \brief Headers of the main subroutines for creating the geometrical structure for the FEM solver.
 *        The subroutines and functions are in the <i>fem_geometry_structure.cpp</i> file.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
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
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "./geometry_structure.hpp"

using namespace std;

/*!
 * \class SortFacesClass
 * \brief Functor, used for a different sorting of the faces than the < operator
 *        of FaceOfElementClass.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class SortFacesClass {
public:
  /*!
   * \brief Constructor of the class. Set the value of nVolElemTot.
   */
   SortFacesClass(unsigned long val_nVolElemTot);

 /*!
  * \brief Destructor of the class. Nothing to be done.
  */
  ~SortFacesClass();

 /*!
  * \brief Operator used for the comparison.
  * \param[in] f0 - First face in the comparison.
  * \param[in] f1 - Second face in the comparison.
  */
  bool operator()(const FaceOfElementClass &f0,
                  const FaceOfElementClass &f1);
private:
  unsigned long nVolElemTot;  /*!< \brief Number of local volume elements . */

  /*!
   * \brief Default constructor of the class. Disabled.
   */
   SortFacesClass(void);

};

/*!
 * \class CPointCompare
 * \brief Helper class used to determine whether two points are identical.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CPointCompare {
public:
  unsigned short nDim;       /*!< \brief Number of spatial dimensions. */
  unsigned long  nodeID;     /*!< \brief The corresponding node ID in the grid. */
  su2double coor[3];         /*!< \brief Coordinates of the point. */
  su2double tolForMatching;  /*!< \brief Tolerance used to determine if points are matching. */

  /*!
   * \brief Constructor of the class. Nothing to be done
   */
  CPointCompare(void);

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  ~CPointCompare(void);

  /*!
   * \brief Copy constructor of the class.
   */
  CPointCompare(const CPointCompare &other);

  /*!
   * \brief Assignment operator of the class.
   */
  CPointCompare& operator=(const CPointCompare &other);

  /*!
   * \brief Less than operator of the class. Needed for the sorting and searching.
   */
  bool operator<(const CPointCompare &other) const;

private:
  /*!
   * \brief Copy function. Needed for the copy constructor and assignment operator.
   */
  void Copy(const CPointCompare &other);
};

/*!
 * \class CVolumeElementFEM
 * \brief Class to store a volume element for the FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CVolumeElementFEM {
public:
  bool elemIsOwned;             /*!< \brief Whether or not this is an owned element. */
  bool JacIsConsideredConstant; /*!< \brief Whether or not the Jacobian of the transformation
                                     to the standard element is considered constant. */

  int rankOriginal;            /*!< \brief The rank where the original volume is stored. For
                                    the original volume, this is simply the current rank. */

  short periodIndexToDonor;    /*!< \brief The index of the periodic transformation to the donor
                                    element. Only for halo elements. A -1 indicates no periodic
                                    transformation. */

  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
  unsigned short nPolyGrid;    /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short nPolySol;     /*!< \brief Polynomial degree for the solution of the element. */
  unsigned short nDOFsGrid;    /*!< \brief Number of DOFs for the geometry of the element. */
  unsigned short nDOFsSol;     /*!< \brief Number of DOFs for the solution of the element. */
  unsigned short nFaces;       /*!< \brief Number of faces of the element. */

  unsigned short indStandardElement; /*!< \brief Index in the vector of standard elements. */

  unsigned long elemIDGlobal;        /*!< \brief Global element ID of this element. */
  unsigned long offsetDOFsSolGlobal; /*!< \brief Global offset of the solution DOFs of this element. */
  unsigned long offsetDOFsSolLocal;  /*!< \brief Local offset of the solution DOFs of this element. */

  vector<bool> JacFacesIsConsideredConstant; /*!< \brief Vector with the booleans whether the Jacobian of the
                                                  transformation to the standard element is constant for the faces. */

  vector<unsigned long> nodeIDsGrid; /*!< \brief Vector with the node IDs of the grid for this element. */

  su2double *metricTerms;            /*!< \brief Pointer to the metric terms in the
                                                 integration points of this element. */

  /*!
   * \brief Constructor of the class. Nothing to be done
   */
  CVolumeElementFEM(void);

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  ~CVolumeElementFEM(void);

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
 * \version 4.1.0 "Cardinal"
 */
class CPointFEM {
public:
  unsigned long globalID;    /*!< \brief The global ID of this point in the grid. */
  short periodIndexToDonor;  /*!< \brief The index of the periodic transformation to the donor
                                  element. Only for halo elements. A -1 indicates no periodic
                                  transformation. */
  su2double coor[3];         /*!< \brief Array with the coordinates of the node. */

  /*!
   * \brief Default constructor of the class. Initialize the coordinates to zero
            to avoid a valgrind warning in two space dimensions.
   */
  CPointFEM(void);

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  ~CPointFEM(void);

  /*!
   * \brief Copy constructor of the class.
   */
  CPointFEM(const CPointFEM &other);

  /*!
   * \brief Assignment operator of the class.
   */
  CPointFEM& operator=(const CPointFEM &other);

  /*!
   * \brief Less than operator of the class. Needed for the sorting.
   */
  bool operator<(const CPointFEM &other) const;

  /*!
   * \brief Equal operator of the class. Needed for the removal of double entities.
   */
  bool operator==(const CPointFEM &other) const;

private:
  /*!
   * \brief Copy function. Needed for the copy constructor and assignment operator.
   */
  void Copy(const CPointFEM &other);
};

/*!
 * \class CInternalFaceElementFEM
 * \brief Class to store an internal face for the FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CInternalFaceElementFEM {
public:
  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */

  unsigned short indStandardElement; /*!< \brief Index in the vector of standard face elements. */

  unsigned long *DOFsGridFaceSide0;   /*!< \brief Pointer to the grid DOFs of side 0 of the face. */
  unsigned long *DOFsGridFaceSide1;   /*!< \brief Pointer to the grid DOFs of side 1 of the face. */
  unsigned long *DOFsSolFaceSide0;    /*!< \brief Pointer to the solution DOFs of side 0 of the face. */
  unsigned long *DOFsSolFaceSide1;    /*!< \brief Pointer to the solution DOFs of side 1 of the face. */

  unsigned long *DOFsGridElementSide0;   /*!< \brief Pointer to the grid DOFs of the element of side 0. */
  unsigned long *DOFsGridElementSide1;   /*!< \brief Pointer to the grid DOFs of the element of side 1. */
  unsigned long *DOFsSolElementSide0;    /*!< \brief Pointer to the solution DOFs of the element of side 0. */
  unsigned long *DOFsSolElementSide1;    /*!< \brief Pointer to the solution DOFs of the element of side 1. */

  /*!
   * \brief Constructor of the class. Nothing to be done.
   */
  CInternalFaceElementFEM(void);

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  ~CInternalFaceElementFEM(void);
};

/*!
 * \class CSurfaceElementFEM
 * \brief Class to store a surface element for the FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CSurfaceElementFEM {
public:
  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
  unsigned short nPolyGrid;    /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short nDOFsGrid;    /*!< \brief Number of DOFs for the geometry of the element. */

  unsigned short indStandardElement; /*!< \brief Index in the vector of standard elements. */

  unsigned long volElemID;         /*!< \brief ID of the corresponding volume element. */
  unsigned long boundElemIDGlobal; /*!< \brief Global ID of this surface element inside
                                        the boundary to which it belongs. */

  vector<unsigned long> nodeIDsGrid; /*!< \brief Vector with the node IDs of the grid for this element.
                                                 In this vector the original sequence of the grid file
                                                 is stored. */

  unsigned long *DOFsGridFace;   /*!< \brief Pointer to the grid DOFs of the face. In principle
                                             the same information as nodeIDsGrid, but the sequence
                                             could be different. */
  unsigned long *DOFsSolFace;    /*!< \brief Pointer to the solution DOFs of the face. */

  unsigned long *DOFsGridElement;   /*!< \brief Pointer to the grid DOFs of the adjacent element. */
  unsigned long *DOFsSolElement;    /*!< \brief Pointer to the solution DOFs of the adjacent element. */

  /*!
   * \brief Constructor of the class. Initialize some variables to avoid a valgrid warning.
   */
  CSurfaceElementFEM(void);

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  ~CSurfaceElementFEM(void);

  /*!
   * \brief Copy constructor of the class.
   */
  CSurfaceElementFEM(const CSurfaceElementFEM &other);

  /*!
   * \brief Assignment operator of the class.
   */
  CSurfaceElementFEM& operator=(const CSurfaceElementFEM &other);

  /*!
   * \brief Less than operator of the class. Needed for the sorting.
   */
  bool operator<(const CSurfaceElementFEM &other) const;

  /*!
   *  \brief Function, which determines a length scale for this surface element.
   *  \return  The relevant length scale for this surface element.
   */
  su2double DetermineLengthScale(vector<CPointFEM> &meshPoints);

  /*!
   *  \brief Function, which determines the corner points of this surface element.
   *  \param[out] nPointsPerFace - Number of corner points of the face.
   *  \param[out] faceConn       - The corner points of the face.
   */
  void GetCornerPointsFace(unsigned short &nPointsPerFace,
                           unsigned long  faceConn[]);

private:
  /*!
   * \brief Copy function. Needed for the copy constructor and assignment operator.
   */
  void Copy(const CSurfaceElementFEM &other);
};

/*!
 * \class CBoundaryFEM
 * \brief Class to store a boundary for the FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CBoundaryFEM {
public:
  string markerTag;  /*!< \brief Marker tag of this boundary. */

  bool periodicBoundary; /*!< \brief Whether or not this boundary is a periodic boundary. */

  vector<CSurfaceElementFEM> surfElem; /*!< \brief Vector of the local surface elements. */

  vector<unsigned long> VecDOFsGridFace; /*!< \brief Storage for the grid DOFs of the faces. */
  vector<unsigned long> VecDOFsSolFace;  /*!< \brief Storage for the solution DOFs of the faces. */

  vector<unsigned long> VecDOFsGridElement; /*!< \brief Storage for the grid DOFs of the adjacent elements. */
  vector<unsigned long> VecDOFsSolElement;  /*!< \brief Storage for the solution DOFs of the adjacent elements. */

  /*!
   * \brief Constructor of the class. Nothing to be done.
   */
  CBoundaryFEM(void);

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  ~CBoundaryFEM(void);
};

/*!
 * \class CMeshFEM
 * \brief Base class for the FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CMeshFEM: public CGeometry {
protected:
  unsigned long nVolElemTot;    /*!< \brief Total number of local volume elements, including halos. */
  unsigned long nVolElemOwned;  /*!< \brief Number of owned local volume elements. */

  vector<CVolumeElementFEM> volElem; /*!< \brief Vector of the local volume elements, including halos. */

  vector<CPointFEM> meshPoints;      /*!< \brief Vector of the points of the FEM mesh. */

  vector<CBoundaryFEM> boundaries;   /*!< \brief Vector of the boundaries of the FEM mesh. */

  vector<unsigned short> rotPerMarkers; /*!< \brief Vector, which contains the indices of the rotational
                                                    periodic markers. */
  vector<vector<unsigned long> > rotPerHalos; /*!< \brief Vector of vector, which contains the indices of
                                                          to halo elements for which a rotationally periodic
                                                          correction must be applied. */

  vector<int> ranksComm;             /*!< \brief Vector of ranks, which this rank exchanges information.
                                                 Self communication is included. */

  vector<vector<unsigned long> > DOFsSend;    /*!< \brief Vector of vector, which contains the DOFs that
                                                          must be sent. Self communication is included. */
  vector<vector<unsigned long> > DOFsReceive; /*!< \brief Vector of vector, which contains the DOFs that
                                                          must be received. Self communication is included. */

  vector<su2double> VecMetricTermsElements;  /*!< \brief Storage for the metric terms of the volume elements. */
  vector<su2double> VecMassMatricesElements; /*!< \brief Storage for the mass matrices of the volume elements. */

public:

  /*!
  * \brief Constructor of the class.
  */
 CMeshFEM(void);

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
  virtual ~CMeshFEM(void);
};

/*!
 * \class CMeshFEM_DG
 * \brief Class which contains all the variables for the DG FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CMeshFEM_DG: public CMeshFEM {
private:
  vector<FEMStandardElementClass> standardElementsSol;  /*!< \brief Vector that contains the standard volume elements
                                                                    used for the solution of the DG solver. */
  vector<FEMStandardElementClass> standardElementsGrid; /*!< \brief Vector that contains the standard volume elements
                                                                    used for the geometry of the DG solver. */

  vector<FEMStandardInternalFaceClass> standardMatchingFacesSol;  /*!< \brief Vector that contains the standard matching
                                                                              internal faces used for the solution of
                                                                              the DG solver. */
  vector<FEMStandardInternalFaceClass> standardMatchingFacesGrid; /*!< \brief Vector that contains the standard matching
                                                                              internal faces used for the geometry of
                                                                              the DG solver. */

  vector<FEMStandardBoundaryFaceClass> standardBoundaryFacesSol;  /*!< \brief Vector that contains the standard boundary
                                                                              faces used for the solution of the DG solver. */
  vector<FEMStandardBoundaryFaceClass> standardBoundaryFacesGrid; /*!< \brief Vector that contains the standard boundary
                                                                              faces used for the geometry of the DG solver. */

  vector<unsigned long> VecDOFsGridFaceSide0; /*!< \brief Storage for the grid DOFs of side 0 of the faces. */
  vector<unsigned long> VecDOFsGridFaceSide1; /*!< \brief Storage for the grid DOFs of side 1 of the faces. */
  vector<unsigned long> VecDOFsSolFaceSide0;  /*!< \brief Storage for the solution DOFs of side 0 of the faces. */
  vector<unsigned long> VecDOFsSolFaceSide1;  /*!< \brief Storage for the solution DOFs of side 1 of the faces. */

  vector<unsigned long> VecDOFsGridElementSide0; /*!< \brief Storage for the grid DOFs of the elements adjacent to side 0. */
  vector<unsigned long> VecDOFsGridElementSide1; /*!< \brief Storage for the grid DOFs of the elements adjacent to side 1. */
  vector<unsigned long> VecDOFsSolElementSide0;  /*!< \brief Storage for the solution DOFs of the elements adjacent to side 0. */
  vector<unsigned long> VecDOFsSolElementSide1;  /*!< \brief Storage for the solution DOFs of the elements adjacent to side 1. */

  vector<CInternalFaceElementFEM> matchingFaces; /*!< \brief Vector of the local matching internal faces. */

public:

 /*!
  * \brief Constructor of the class.
  */
 CMeshFEM_DG(void);

  /*!
  * \overload
  * \brief Redistributes the grid over the ranks and creates the halo layer.
  * \param[in] geometry - The linear distributed grid that must be redistributed.
  * \param[in] config   - Definition of the particular problem.
  */
  CMeshFEM_DG(CGeometry *geometry, CConfig *config);

 /*!
  * \brief Destructor of the class.
  */
  ~CMeshFEM_DG(void);

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
  * \brief Function, which computes the metric terms of the
           surface elements.
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
  void SetSendReceive(CConfig *config);

private:
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
};

#include "fem_geometry_structure.inl"
