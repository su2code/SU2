/*!
 * \file CFEMStandardQuadGrid.hpp
 * \brief Class for the FEM quadrilateral standard element for the grid.
 *        The functions are in the <i>CFEMStandardQuadGrid.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
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

#include "CFEMStandardQuad.hpp"

/*!
 * \class CFEMStandardQuadGrid
 * \brief Class which defines the variables and methods for the
 *        quadrilateral standard element for the grid.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
 */
class CFEMStandardQuadGrid final: public CFEMStandardQuad {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardQuadGrid() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly       - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact  - Polynomial order that must be integrated exactly
   *                              by the integration rule.
   * \param[in] val_surfElement - True if this element is a surface element,
   *                              False if this element is a volume element (2D simulation)
   */
  CFEMStandardQuadGrid(const unsigned short val_nPoly,
                       const unsigned short val_orderExact,
                       const bool           val_surfElement);

  /*!
   * \overload
   * \param[in] val_nPolyGrid   - Polynomial degree of the grid for this element.
   * \param[in] val_nPolyGrid   - Polynomial degree of the solution for this element.
   * \param[in] val_orderExact  - Polynomial order that must be integrated exactly
   *                              by the integration rule.
   * \param[in] val_locGridDOFs - Location of the grid DOFS, either LGL or equidistant.
   */
  CFEMStandardQuadGrid(const unsigned short val_nPolyGrid,
                       const unsigned short val_nPolySol,
                       const unsigned short val_orderExact,
                       const unsigned short val_locGridDOFs);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  ~CFEMStandardQuadGrid() = default;

  /*!
   * \brief Function, which computes the coordinates in the integration points.
   * \param[in]  LGLDistribution - Whether or not the LGL node distribution must be used.
   * \param[in]  matCoorDOF - Matrix that contains the coordinates of the grid DOFs.
   * \param[out] matCoorInt - Matrix that contains the coordinates of the integration
   *                          points.
   */
  void CoorIntPoints(const bool                LGLDistribution,
                     ColMajorMatrix<su2double> &matCoorDOF,
                     ColMajorMatrix<su2double> &matCoorInt) override;

  /*!
   * \brief Function, which computes the derivatives of the coordinates in the
   *        integration points.
   * \param[in]  LGLDistribution - Whether or not the LGL node distribution must be used.
   * \param[in]  matCoor         - Matrix that contains the coordinates of the grid DOFs.
   * \param[out] matDerCoor      - Vector of matrices to store the derivatives of the coordinates.
   */
  void DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                ColMajorMatrix<su2double>          &matCoor,
                                vector<ColMajorMatrix<su2double> > &matDerCoor) override;

  /*!
   * \brief Function, which computes the derivatives of the coordinates in the
   *        solution DOFs.
   * \param[in]  matCoor    - Matrix that contains the coordinates of the grid DOFs.
   * \param[out] matDerCoor - Vector of matrices to store the derivatives of the coordinates.
   */
  void DerivativesCoorSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                              vector<ColMajorMatrix<su2double> > &matDerCoor) override;

  /*!
   * \brief Function, that returns the number of different face types
   *        occuring in this volume element.
   * \return The number of different face types of the volume element.
   */
  unsigned short GetnFaceTypes(void) const override {return 1;}

  /*!
   * \brief Function, that returns the number of solution DOFs.
   * \return The number of solution DOFs of the volume element.
   */
  unsigned short GetNSolDOFs(void) const override {
    const unsigned short nSol1D = rLineSolDOFs.size();
    return nSol1D*nSol1D;
  }

  /*!
   * \brief Function that returns the VTK type for the given face type index.
   * \param[in] ind - Index of the face type for which the VTK type must be returned.
   * \return The VTK type of the given face type.
   */
  unsigned short GetVTK_TypeFace(unsigned short ind) const override {return LINE;}

  /*!
   * \brief Function, which estimates the amount of work for a boundary quad face.
   *        This information is used to determine a well balanced partition.
   * \param[in] config   - Object, which contains the input parameters.
   * \param[in] elemType - Type of the volume element adjacent to this boundary face.
   * \return The work estimate for a boundary face of this type.
   */
  passivedouble WorkEstimateBoundaryFace(CConfig              *config,
                                         const unsigned short elemType) override;

  /*!
   * \brief Function, which estimates the amount of work for an internal quad face.
   *        This information is used to determine a well balanced partition.
   * \param[in] config    - Object, which contains the input parameters.
   * \param[in] elemType0 - Type of the volume element adjacent to side 0 of this face.
   * \param[in] nPoly0    - Polynomial degree used in elemType0.
   * \param[in] elemType1 - Type of the volume element adjacent to side 1 of this face.
   * \param[in] nPoly1    - Polynomial degree used in elemType1.
   * \return The work estimate for an internal face of this type.
   */
  passivedouble WorkEstimateInternalFace(CConfig              *config,
                                         const unsigned short elemType0,
                                         const unsigned short nPoly0,
                                         const unsigned short elemType1,
                                         const unsigned short nPoly1) override;

  /*!
   * \brief Function, which estimates the amount of work for a volume quadrilateral.
   *        This information is used to determine a well balanced partition.
   *        The work of the surface integral in DG is not included.
   * \param[in] config - Object, which contains the input parameters.
   * \return The work estimate for the volume for this type of element.
   */
  passivedouble WorkEstimateVolume(CConfig *config) override;

  /*!
   * \brief Function, which estimates the amount of work for a boundary quad
   *        when wall functions are applied.
   * \param[in] config    - Object, which contains the input parameters.
   * \param[in] nPointsWF - Number of points in wall function treatment.
   * \param[in] elemType  - Type of the volume element adjacent to this boundary face.
   * \return The work estimate for a boundary face of this type.
   */
  passivedouble WorkEstimateWallFunctions(CConfig              *config,
                                          const unsigned short nPointsWF,
                                          const unsigned short elemType) override;
private:

  unsigned short nDim; /*!< \brief Number of space dimensions. For nDim = 2 this standard element is
                                   a volume element, while for nDim = 3 it is a surface element. */

  vector<passivedouble> rLineSolDOFs; /*!< \brief 1D parametric coordinates of the nodal solution DOFs. */

  ColMajorMatrix<passivedouble> lagBasisLineIntEqui; /*!< \brief The values of the 1D Lagrangian basis functions
                                                                 in the integration points for the equidistant
                                                                 point distribution. */
  ColMajorMatrix<passivedouble> lagBasisLineIntLGL;  /*!< \brief The values of the 1D Lagrangian basis functions
                                                                 in the integration points for the LGL
                                                                 point distribution. */

  ColMajorMatrix<passivedouble> derLagBasisLineIntEqui; /*!< \brief The values of the derivatives of the 1D Lagrangian
                                                                    basis functions in the integration points for the
                                                                    equidistant point distribution. */
  ColMajorMatrix<passivedouble> derLagBasisLineIntLGL;  /*!< \brief The values of the derivatives of the 1D Lagrangian
                                                                    basis functions in the integration points for the
                                                                    LGL point distribution. */

  ColMajorMatrix<passivedouble> lagBasisLineSolDOFs;    /*!< \brief The values of the 1D Lagrangian basis functions
                                                                    in the nodal solution DOFs. */
  ColMajorMatrix<passivedouble> derLagBasisLineSolDOFs; /*!< \brief The values of the derivatives of the 1D Lagrangian
                                                                    basis functions in the nodal solution DOFs. */

  ColMajorMatrix<passivedouble> hesLagBasisLineInt;     /*!< \brief The values of the 2nd derivatives, Hessian, of the
                                                                    Lagrangian basis function in the integration points
                                                                    for the point distribution that is in use. */

  TPI2D TensorProductDataVolSolDOFs = nullptr; /*!< \brief Function pointer to carry out the tensor product
                                                           to compute the data in the nodal solution DOFs. */

  /*!
   * \brief Function, which creates the local grid connectivities of the faces
   *        of the volume element.
   */
  void LocalGridConnFaces(void);

  /*!
   * \brief Function, which creates the connectivity of the linear sub-elements when the
   *        high order element is split in such elements.
   */
  void SubConnLinearElements(void);

  /*!
   * \brief Function, which serves as an interface to carry out the tensor product C = A*B
   *        to obtain the data in the solution DOFs.
   * \param[in]  N      - Number of values to be computed in the integration points.
   * \param[in]  Ai     - 1D matrix for the i-component of the A tensor.
   * \param[in]  Aj     - 1D matrix for the j-component of the A tensor.
   * \param[in]  B      - B tensor stored as a matrix.
   * \param[out] C      - C tensor stored as a matrix.
   * \param[out] config - Object used for the timing of the tensor product call.
   */
  void TensorProductSolDOFs(const int                           N,
                            const ColMajorMatrix<passivedouble> &Ai,
                            const ColMajorMatrix<passivedouble> &Aj,
                            const ColMajorMatrix<su2double>     &B,
                            ColMajorMatrix<su2double>           &C,
                            const CConfig                       *config);
};
