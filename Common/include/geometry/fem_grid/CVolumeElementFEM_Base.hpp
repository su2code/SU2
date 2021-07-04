/*!
 * \file CVolumeElementFEM_Base.hpp
 * \brief Base class for a volume element for the FEM solver.
 *        The implementations are in the <i>CVolumeElementFEM_Base.cpp</i> file.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
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

#include "../../fem/CFEMStandardElementBase.hpp"
#include "../../adt/CADTElemClass.hpp"

using namespace std;

/*--- Forward declaration. ---*/
class CPointFEM;

/*!
 * \class CVolumeElementFEM_Base
 * \brief Base class to store a volume element for the FEM solver.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
 */
class CVolumeElementFEM_Base {
public:
  bool elemIsOwned;             /*!< \brief Whether or not this is an owned element. */
  bool JacIsConsideredConstant; /*!< \brief Whether or not the Jacobian of the transformation
                                            to the standard element is considered constant. */

  int rankOriginal;            /*!< \brief The rank where the original volume is stored. For
                                           the owned volumes, this is simply the current rank. */

  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
  unsigned short nPolyGrid;    /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short nPolySol;     /*!< \brief Polynomial degree for the solution of the element. */
  unsigned short nDOFsGrid;    /*!< \brief Number of DOFs for the geometry of the element. */
  unsigned short nDOFsSol;     /*!< \brief Number of DOFs for the solution of the element. */
  unsigned short nFaces;       /*!< \brief Number of faces of the element. */

  unsigned long elemIDGlobal;        /*!< \brief Global element ID of this element. */

  vector<bool> JacFacesIsConsideredConstant; /*!< \brief Vector with the booleans whether the Jacobian of the
                                                         transformation to the standard element is constant
                                                         for the faces. */
  vector<bool> ElementOwnsFaces;             /*!< \brief Vector with the booleans whether the element is the
                                                         owner of the faces. */

  vector<unsigned long> nodeIDsGrid; /*!< \brief Vector with the node IDs of the grid for this element. */

  su2double lenScale;                /*!< \brief Length scale of the element. */
  su2double volume;                  /*!< \brief Volume of the element. */
  su2double avgJacobian;             /*!< \brief Average value of the Jacobian of the element. */

  ColMajorMatrix<su2double> coorGridDOFs;         /*!< \brief Coordinates of the grid DOFs of this element. */

  su2activevector JacobiansInt;      /*!< \brief The Jacobians in the integration points of this element. */ 
  su2activevector JacobiansSolDOFs;  /*!< \brief The Jacobians in the nodal solution DOFs of this element. */

  vector<ColMajorMatrix<su2double> > metricTermsInt;     /*!< \brief The metric terms in the integration
                                                                     points of this element. */
  vector<ColMajorMatrix<su2double> > metricTermsSolDOFs; /*!< \brief The metric terms in the nodal
                                                                     solution DOFs of this element. */

  vector<ColMajorMatrix<su2double> > metricTerms2ndDerInt; /*!< \brief The metric terms needed for the computation
                                                                       of the 2nd derivatives in the integration
                                                                       points. Only determined when needed (ADER-DG
                                                                       with non-aliased predictor for the 
                                                                       Navier-Stokes equations). */

  ColMajorMatrix<su2double> coorIntegrationPoints;  /*!< \brief Coordinates of the integration points of this element. */
  ColMajorMatrix<su2double> gridVelocitiesInt;      /*!< \brief Grid velocities of the integration points of this element. */
  su2activevector           wallDistanceInt;        /*!< \brief The wall distance to the viscous walls for
                                                                the integration points of this element. */

  ColMajorMatrix<su2double> coorSolDOFs;            /*!< \brief Coordinates of the nodal solution DOFs of this element.
                                                                It only differs from coorGridDOFs if the polynomial degree
                                                                of the grid and solution differ. */
  ColMajorMatrix<su2double> gridVelocitiesSolDOFs;  /*!< \brief Grid velocities of the nodal solution DOFs of this element. */
  su2activevector           wallDistanceSolDOFs;    /*!< \brief The wall distance to the viscous walls for
                                                                the nodal solution DOFs of this element. */    

  CFEMStandardElementBase *standardElemGrid = nullptr; /*!< \brief Pointer to the standard element for the grid. */

  /*!
   * \brief Function, which computes the wall distances.
   * \param[in] WallADT - The ADT to compute the wall distances.
   * \param[in] nDim    - Number of spatial dimensions.
   */
  void ComputeWallDistance(CADTElemClass        *WallADT,
                           const unsigned short nDim);

  /*!
   * \brief Function, which computes the derivatives of the metric terms
   *         in the integration points.
   * \param[in] nDim - Number of spatial dimensions.
   */
  void DerMetricTermsIntegrationPoints(const unsigned short nDim);

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
                               unsigned long  faceConn[6][4]) const;

  /*!
   * \brief Function, which initializes the grid velocities of this volume element.
   * \param[in] nDim - Number of spatial dimensions.
   */
  void InitGridVelocities(const unsigned short nDim);

  /*!
   * \brief Function, which computes the metric terms in the integration points.
   * \param[in] LGLDistribution - Whether or not the LGL node distribution must be used.
   * \param[in] nDim            - Number of spatial dimensions.
   * \return  True of all Jacobians are positive and negative otherwise.
   */
  bool MetricTermsIntegrationPoints(const bool           LGLDistribution,
                                    const unsigned short nDim);

  /*!
   * \brief Function, which computes the metric terms in the nodal solution DOFs.
   * \param[in] nDim - Number of spatial dimensions.
   * \return  True of all Jacobians are positive and negative otherwise.
   */
  bool MetricTermsSolDOFs(const unsigned short nDim);

  /*!
   * \brief Function, which sets the coordinates of the grid DOFs.
   * \param[in] nDim       - Number of spatial dimensions.
   * \param[in] meshPoints - Vector, which contains the coordinates
   *                         of all points in the current partition.
   */
  void SetCoorGridDOFs(const unsigned short    nDim,
                       const vector<CPointFEM> &meshPoints);

  /*!
   * \brief Function, which sets the wall distances to the given value.
   * \param[in] val - Value to which the wall distance must be set.
   */
  void SetWallDistance(su2double val);
};
