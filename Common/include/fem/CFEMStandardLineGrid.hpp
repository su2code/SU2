/*!
 * \file CFEMStandardLineGrid.hpp
 * \brief Class for the FEM line standard element for the grid.
 *        The functions are in the <i>CFEMStandardaceLineGrid.cpp</i> file.
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

#include "CFEMStandardLine.hpp"

/*!
 * \class CFEMStandardLineGrid
 * \brief Class which defines the variables and methods for the
 *        line standard element for the grid.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 */
class CFEMStandardLineGrid final: public CFEMStandardLine {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardLineGrid() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardLineGrid(const unsigned short val_nPoly,
                       const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  ~CFEMStandardLineGrid() = default;

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
   * \brief Function, which estimates the amount of work for a boundary line face.
   *        This information is used to determine a well balanced partition.
   * \param[in] config   - Object, which contains the input parameters.
   * \param[in] elemType - Type of the volume element adjacent to this boundary face.
   * \return The work estimate for a boundary face of this type.
   */
  passivedouble WorkEstimateBoundaryFace(CConfig              *config,
                                         const unsigned short elemType) override;

  /*!
   * \brief Function, which estimates the amount of work for an internal line face.
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
   * \brief Function, which estimates the amount of work for a boundary line
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

  /*!
   * \brief Function, which creates the connectivity of the linear sub-elements when the
   *        high order element is split in such elements.
   */
  void SubConnLinearElements(void);
};
