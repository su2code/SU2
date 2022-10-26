/*!
 * \file CFEMStandardTriPartition.hpp
 * \brief Class for the FEM triangle standard element
 *        used during the partitioning.
 *        The functions are in the <i>CFEMStandardTriPartition.cpp</i> file.
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

#include "CFEMStandardTriBase.hpp"

/*!
 * \class CFEMStandardTriPartition
 * \brief Class which defines the variables and methods for the
 *        triangle standard element used in the partitioning
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
 */
class CFEMStandardTriPartition final: public CFEMStandardTriBase {

public:
  /*-----------------------------------------------------------------------------------*/
  /*---                     Constructors and destructors.                           ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardTriPartition() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly       - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact  - Polynomial order that must be integrated exactly
   *                              by the integration rule.
   * \param[in] val_surfElement - True if this element is a surface element,
   *                              False if this element is a volume element (2D simulation)
   */
  CFEMStandardTriPartition(const unsigned short val_nPoly,
                           const unsigned short val_orderExact,
                           const bool           val_surfElement);

  /*!
   * \brief Destructor.
   */
  ~CFEMStandardTriPartition();

  /*-----------------------------------------------------------------------------------*/
  /*---                  Inline public member functions.                            ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, that returns the number of different face types
   *        occuring in this volume element.
   * \return The number of different face types of the volume element.
   */
  inline unsigned short GetnFaceTypes(void) const override {return 1;}

  /*!
   * \brief Function that returns the VTK type for the given face type index.
   * \param[in] ind - Index of the face type for which the VTK type must be returned.
   * \return The VTK type of the given face type.
   */
  inline unsigned short GetVTK_TypeFace(unsigned short ind) const override {return LINE;}

  /*-----------------------------------------------------------------------------------*/
  /*---                       Public member functions.                              ---*/
  /*-----------------------------------------------------------------------------------*/

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
   * \brief Function, which estimates the amount of work for a boundary triangular face.
   *        This information is used to determine a well balanced partition.
   * \param[in] config   - Object, which contains the input parameters.
   * \param[in] elemType - Type of the volume element adjacent to this boundary face.
   * \return The work estimate for a boundary face of this type.
   */
  passivedouble WorkEstimateBoundaryFace(CConfig              *config,
                                         const unsigned short elemType) override;

  /*!
   * \brief Function, which estimates the amount of work for an internal triangular face.
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
   * \brief Function, which estimates the amount of work for a volume triangle.
   *        This information is used to determine a well balanced partition.
   *        The work of the surface integral in DG is not included.
   * \param[in] config - Object, which contains the input parameters.
   * \return The work estimate for the volume for this type of element.
   */
  passivedouble WorkEstimateVolume(CConfig *config) override;

  /*!
   * \brief Function, which estimates the amount of work for a boundary triangle
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

  vector<passivedouble> rTriangleDOFsEqui; /*!< \brief Parametric r-coordinates of the triangle grid
                                                       DOFs when equidistant spacing is used. */
  vector<passivedouble> sTriangleDOFsEqui; /*!< \brief Parametric s-coordinates of the triangle grid
                                                       DOFs when equidistant spacing is used. */
  vector<passivedouble> rTriangleDOFsLGL;  /*!< \brief Parametric r-coordinates of the triangle grid
                                                       DOFs when the LGL distribution is used. */
  vector<passivedouble> sTriangleDOFsLGL;  /*!< \brief Parametric s-coordinates of the triangle grid
                                                       DOFs when the LGL distribution is used. */

  ColMajorMatrix<passivedouble> lagBasisIntEqui; /*!< \brief The values of the Lagrangian basis functions
                                                             in the integration points for the equidistant
                                                             point distribution. */
  ColMajorMatrix<passivedouble> lagBasisIntLGL;  /*!< \brief The values of the Lagrangian basis functions
                                                             in the integration points for the LGL
                                                             point distribution. */

  vector<ColMajorMatrix<passivedouble> > derLagBasisIntEqui; /*!< \brief The values of the derivatives of the Lagrangian
                                                                         basis functions in the integration points for the
                                                                         equidistant point distribution. It is a vector,
                                                                         because there are derivatives in two directions. */
  vector<ColMajorMatrix<passivedouble> > derLagBasisIntLGL;  /*!< \brief The values of the derivatives of the Lagrangian
                                                                         basis functions in the integration points for the
                                                                         LGL point distribution. It is a vector, because
                                                                         there are derivatives in two directions. */

  void *jitterDOFs2Int = nullptr;      /*!< \brief Pointer to the data for the jitted gemm function
                                                   to compute data in the integration points. */
  dgemm_jit_kernel_t gemmDOFs2Int;     /*!< \brief Pointer to the function to carry out jitterDOFs2Int. */
};
