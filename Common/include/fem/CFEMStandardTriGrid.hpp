/*!
 * \file CFEMStandardTriGrid.hpp
 * \brief Class for the FEM triangle standard element for the grid.
 *        The functions are in the <i>CFEMStandardTriGrid.cpp</i> file.
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

#include "CFEMStandardTri.hpp"

/*!
 * \class CFEMStandardTriGrid
 * \brief Class which defines the variables and methods for the
 *        triangle standard element for the grid.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
class CFEMStandardTriGrid final: public CFEMStandardTri {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardTriGrid() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly       - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact  - Polynomial order that must be integrated exactly
   *                              by the integration rule.
   * \param[in] val_surfElement - True if this element is a surface element,
   *                              False if this element is a volume element (2D simulation)
   */
  CFEMStandardTriGrid(const unsigned short val_nPoly,
                      const unsigned short val_orderExact,
                      const bool           val_surfElement);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  ~CFEMStandardTriGrid() = default;

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
   * \brief Function, that returns the number of different face types
   *        occuring in this volume element.
   * \return The number of different face types of the volume element.
   */
  unsigned short GetnFaceTypes(void) const override {return 1;}

  /*!
   * \brief Function that returns the VTK type for the given face type index.
   * \param[in] ind - Index of the face type for which the VTK type must be returned.
   * \return The VTK type of the given face type.
   */
  unsigned short GetVTK_TypeFace(unsigned short ind) const override {return LINE;}

private:

  unsigned short nDim; /*!< \brief Number of space dimensions. For nDim = 2 this standard element is
                                   a volume element, while for nDim = 3 it is a surface element. */

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

  /*!
   * \brief Function, which computes the values of the derivatives of the Lagrangian
   *        basis functions of a triangle in the integration points for the given
   *        location of the DOFs.
   * \param[in]  rDOFs  - Vector, which contains the parametric r-locations of the DOFs.
   * \param[in]  sDOFs  - Vector, which contains the parametric s-locations of the DOFs.
   * \param[out] derLag - Matrix, which contains the values of derivatives of all the
   *                      Lagrangian basis functions in all the integration points.
   */
  void DerLagBasisIntPointsTriangle(const vector<passivedouble>            &rDOFs,
                                    const vector<passivedouble>            &sDOFs,
                                    vector<ColMajorMatrix<passivedouble> > &derLag);

  /*!
   * \brief Function, which computes the values of the Lagrangian basis functions
   *        of a triangle in the integration points for the given location of the DOFs.
   * \param[in]  rDOFs - Vector, which contains the parametric r-locations of the DOFs.
   * \param[in]  sDOFs - Vector, which contains the parametric s-locations of the DOFs.
   * \param[out] lag   - Matrix, which contains the values of all the Lagrangian
   *                     basis functions in all the integration points.
   */
  void LagBasisIntPointsTriangle(const vector<passivedouble>   &rDOFs,
                                 const vector<passivedouble>   &sDOFs,
                                 ColMajorMatrix<passivedouble> &lag);

  /*!
   * \brief Function, which computes the gradients of the Vandermonde matrix for a standard triangle.
   * \param[in]  r   - Parametric r-coordinates for which the Vandermonde matrix must be computed.
   * \param[in]  s   - Parametric s-coordinates for which the Vandermonde matrix must be computed
   * \param[out] VDr - Matrix to store the derivative in r-direction of the Vandermonde matrix
   *                   in all r,s-locations.
   * \param[out] VDs - Matrix to store the derivative in s-direction of the Vandermonde matrix
   *                   in all r,s-locations.
   */
  void GradVandermondeTriangle(const vector<passivedouble>   &r,
                               const vector<passivedouble>   &s,
                               ColMajorMatrix<passivedouble> &VDr,
                               ColMajorMatrix<passivedouble> &VDs);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard triangle.
   * \param[in]  r  - Parametric r-coordinates for which the Vandermonde matrix must be computed.
   * \param[in]  s  - Parametric s-coordinates for which the Vandermonde matrix must be computed
   * \param[out] V  - Matrix to store the Vandermonde matrix in all r,s-locations.
   */
  void VandermondeTriangle(const vector<passivedouble>   &r,
                           const vector<passivedouble>   &s,
                           ColMajorMatrix<passivedouble> &V);
};
