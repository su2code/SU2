/*!
 * \file CSurfaceElementFEM.hpp
 * \brief Class definition for a surface element for the FEM solver.
 *        The implementations are in the <i>CSurfaceElementFEM.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
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

using namespace std;

/*--- Forward declaration. ---*/
class CVolumeElementFEM_DG;

/*!
 * \class CSurfaceeSurfaceFEM
 * \brief Class to store a surface element for the FEM solver.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 */
class CSurfaceElementFEM final {
public:
  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
  unsigned short nPolyGrid;    /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short nDOFsGrid;    /*!< \brief Number of DOFs for the geometry of the element. */

  unsigned long volElemID;         /*!< \brief ID of the corresponding volume element. */
  unsigned long boundElemIDGlobal; /*!< \brief Global ID of this surface element inside
                                               the boundary to which it belongs. */

  vector<unsigned long> nodeIDsGrid; /*!< \brief Vector with the node IDs of the grid for this element. */

  su2activevector JacobiansFace;               /*!< \brief The Jacobians in the integration points of the face. */
  ColMajorMatrix<su2double> metricNormalsFace; /*!< \brief The normals in the integration points of the face.
                                                           The normals point from side 0 out of the element. */

  vector<ColMajorMatrix<su2double> > metricCoorDerivFace;  /*!< \brief The terms drdx, dsdx, etc. of side 0 in the
                                                                       integration points of the face. */

  ColMajorMatrix<su2double> coorIntegrationPoints;  /*!< \brief Coordinates of the integration points of this face. */
  ColMajorMatrix<su2double> gridVelocities;         /*!< \brief Grid velocities of the integration points of this face. */

  su2activevector wallDistance;     /*!< \brief The wall distance to the viscous walls for
                                                the integration points of this face. */

  vector<unsigned long>  donorsWallFunction;                   /*!< \brief Local element IDs of the donors for the wall
                                                                           function treatment. These donors can be halo's. */
  vector<unsigned short> nIntPerWallFunctionDonor;             /*!< \brief The number of integration points per donor
                                                                           element for the wall function treatment. */
  vector<unsigned short> intPerWallFunctionDonor;              /*!< \brief The integration points per donor element
                                                                           for the wall function treatment. */
  vector<ColMajorMatrix<passivedouble> > matWallFunctionDonor; /*!< \brief Matrices, which store the interpolation coefficients
                                                                           for the donors of the integration points.*/

  CFEMStandardElementBase *standardElemGrid = nullptr; /*!< \brief Pointer to the standard element for the grid. */
  CFEMStandardElementBase *standardElemFlow = nullptr; /*!< \brief Pointer to the standard element for the
                                                                   standard flow solution variables. */
  CFEMStandardElementBase *standardElemP    = nullptr; /*!< \brief Pointer to the standard element for the
                                                                   pressure for an incompressible flow. */

  /*!
   * \brief Less than operator of the class. Needed for the sorting.
            The criterion for comparison is the corresponding (local) volume ID.
   */
  bool operator<(const CSurfaceElementFEM &other) const { return volElemID < other.volElemID; }

  /*!
   *  \brief Function, which determines the corner points of this surface element.
   *  \param[out] nPointsPerFace - Number of corner points of the face.
   *  \param[out] faceConn       - The corner points of the face.
   */
  void GetCornerPointsFace(unsigned short &nPointsPerFace,
                           unsigned long  faceConn[]);

  /*!
   * \brief Function, which initializes the grid velocities of this face.
   * \param[in] nDim - Number of spatial dimensions.
   */
  void InitGridVelocities(const unsigned short nDim);

  /*!
   * \brief Function, which computes the metric terms in the integration points.
   * \param[in] nDim    - Number of spatial dimensions.
   * \param[in] volElem - The volume elements of the grid.
   */
  void MetricTermsIntegrationPoints(const unsigned short         nDim,
                                    vector<CVolumeElementFEM_DG> &volElem);
};
