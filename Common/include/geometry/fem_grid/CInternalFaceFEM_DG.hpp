/*!
 * \file CInternalFaceFEM_DG.hpp
 * \brief Class for an internal face element for the DG-FEM solver.
 *        The implementations are in the <i>CInternalFaceFEM_DG.cpp</i> file.
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

#include "../../fem/CFEMStandardInternalFaceGrid.hpp"
#include "../../fem/CFEMStandardInternalFaceSol.hpp"
#include "../../adt/CADTElemClass.hpp"

using namespace std;

/*--- Forward declaration. ---*/
class CVolumeElementFEM_DG;

/*!
 * \class CInternalFaceElementFEM
 * \brief Class to store an internal face for the FEM solver.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
 */
class CInternalFaceFEM_DG final {
public:
  unsigned long elemID0;         /*!< \brief Element ID adjacent to side 0 of the face. */
  unsigned long elemID1;         /*!< \brief Element ID adjacent to side 1 of the face. */

  su2activevector JacobiansFace;               /*!< \brief The Jacobians in the integration points of the face. */
  ColMajorMatrix<su2double> metricNormalsFace; /*!< \brief The normals in the integration points of the face.
                                                           The normals point from side 0 to side 1. */

  vector<ColMajorMatrix<su2double> > metricCoorDerivFace0;  /*!< \brief The terms drdx, dsdx, etc. of side 0 in the
                                                                        integration points of the face. */
  vector<ColMajorMatrix<su2double> > metricCoorDerivFace1;  /*!< \brief The terms drdx, dsdx, etc. of side 1 in the
                                                                        integration points of the face. */

  ColMajorMatrix<su2double> coorIntegrationPoints;  /*!< \brief Coordinates of the integration points of this face. */
  ColMajorMatrix<su2double> gridVelocities;         /*!< \brief Grid velocities of the integration points of this face. */

  su2activevector wallDistance;      /*!< \brief The wall distance to the viscous walls for
                                                 the integration points of this face. */

  ColMajorMatrix<su2double> resDOFsSide0;         /*!< \brief The contribution of this face to the residual
                                                              in the solution DOFs of the element on side0. */
  ColMajorMatrix<su2double> resDOFsSide1;         /*!< \brief The contribution of this face to the residual
                                                              in the solution DOFs of the element on side1. */

  CFEMStandardInternalFaceGrid *standardElemGrid = nullptr; /*!< \brief Pointer to the standard element for the grid. */
  CFEMStandardInternalFaceSol  *standardElemFlow = nullptr; /*!< \brief Pointer to the standard element for the
                                                                        standard flow solution variables. */
  CFEMStandardInternalFaceSol  *standardElemP    = nullptr; /*!< \brief Pointer to the standard element for the
                                                                        pressure for an incompressible flow. */

  /*!
   * \brief Function, which allocate the memory for the residuals.
   * \param[in] config - Definition of the particular problem.
   * \param[in] nVar   - Number of flow variables.
   */
  void AllocateResiduals(CConfig        *config,
                         unsigned short nVar);

  /*!
   * \brief Function, which computes the gradients w.r.t. the parametric coordinates
   *        of the solution on side 0 in the integration points.
   * \param[in] volElem - Pointer to the local volume elements.
   * \return  A reference to the gradients of the solution on side 0
   *          in the integration points.
   */
  vector<ColMajorMatrix<su2double> > &ComputeGradSolSide0IntPoints(CVolumeElementFEM_DG *volElem);

  /*!
   * \brief Function, which computes the gradients w.r.t. the parametric coordinates
   *        of the solution on side 1 in the integration points.
   * \param[in] volElem - Pointer to the local volume elements.
   * \return  A reference to the gradients of the solution on side 1
   *          in the integration points.
   */
  vector<ColMajorMatrix<su2double> > &ComputeGradSolSide1IntPoints(CVolumeElementFEM_DG *volElem);

  /*!
   * \brief Function, which computes the solution on side 0 in the integration points.
   * \param[in] volElem - Pointer to the local volume elements.
   * \return  A reference to the solution on side 0 in the integration points.
   */
  ColMajorMatrix<su2double> &ComputeSolSide0IntPoints(CVolumeElementFEM_DG *volElem);

  /*!
   * \brief Function, which computes the solution on side 1 in the integration points.
   * \param[in] volElem - Pointer to the local volume elements.
   * \return  A reference to the solution on side 1 in the integration points.
   */
  ColMajorMatrix<su2double> &ComputeSolSide1IntPoints(CVolumeElementFEM_DG *volElem);

  /*!
   * \brief Function, which computes the wall distances.
   * \param[in] WallADT - The ADT to compute the wall distances.
   * \param[in] nDim    - Number of spatial dimensions.
   */
  void ComputeWallDistance(CADTElemClass        *WallADT,
                           const unsigned short nDim);

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

  /*!
   * \brief Function, which adds to the residual of element on side 0
   *        the contribution coming from the multiplication with the
   *        basis functions.
   * \param[in] scalarDataInt - The scalar data in the integration points
   *                            to be multiplied by the basis functions.
   */
  void ResidualBasisFunctionsSide0(ColMajorMatrix<su2double> &scalarDataInt);

  /*!
   * \brief Function, which adds to the residual of element on side 1
   *        the contribution coming from the multiplication with the
   *        basis functions.
   * \param[in] scalarDataInt - The scalar data in the integration points
   *                            to be multiplied by the basis functions.
   */
  void ResidualBasisFunctionsSide1(ColMajorMatrix<su2double> &scalarDataInt);

  /*!
   * \brief Function, which adds to the residual of element on side 0 the contribution
   *        coming from the multiplication with the divergence of the basis functions.
   * \param[in] vectorDataInt - The vector data in the integration points to be
   *                            multiplied by the gradient of the basis functions.
   */
  void ResidualGradientBasisFunctionsSide0(vector<ColMajorMatrix<su2double> > &vectorDataInt);

  /*!
   * \brief Function, which adds to the residual of element on side 1 the contribution
   *        coming from the multiplication with the divergence of the basis functions.
   * \param[in] vectorDataInt - The vector data in the integration points to be
   *                            multiplied by the gradient of the basis functions.
   */
  void ResidualGradientBasisFunctionsSide1(vector<ColMajorMatrix<su2double> > &vectorDataInt);

  /*!
   * \brief Function, which sets the wall distances to the given value.
   * \param[in] val - Value to which the wall distance must be set.
   */
  void SetWallDistance(su2double val);
};
