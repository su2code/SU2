/*!
 * \file CFEMStandardPrismGrid.hpp
 * \brief Class for the FEM prism standard element for the grid.
 *        The functions are in the <i>CFEMStandardPrismGrid.cpp</i> file.
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

#include "CFEMStandardPrism.hpp"

/*!
 * \class CFEMStandardPrismGrid
 * \brief Class which defines the variables and methods for the
 *        prism standard element for the grid.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
 */
class CFEMStandardPrismGrid final: public CFEMStandardPrism {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardPrismGrid() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardPrismGrid(const unsigned short val_nPoly,
                        const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  ~CFEMStandardPrismGrid() = default;

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
  unsigned short GetnFaceTypes(void) const override {return 2;}

  /*!
   * \brief Function that returns the VTK type for the given face type index.
   * \param[in] ind - Index of the face type for which the VTK type must be returned.
   * \return The VTK type of the given face type.
   */
  unsigned short GetVTK_TypeFace(unsigned short ind) const override {
    if(ind == 0) return TRIANGLE;
    else         return QUADRILATERAL;
  }

private:

  ColMajorMatrix<passivedouble> lagBasisIntEqui; /*!< \brief The values of the Lagrangian basis functions
                                                             in the integration points for the equidistant
                                                             point distribution. */
  ColMajorMatrix<passivedouble> lagBasisIntLGL;  /*!< \brief The values of the Lagrangian basis functions
                                                             in the integration points for the LGL
                                                             point distribution. */

  vector<ColMajorMatrix<passivedouble> > derLagBasisIntEqui; /*!< \brief The values of the derivatives of the Lagrangian
                                                                         basis functions in the integration points for the
                                                                         equidistant point distribution. It is a vector,
                                                                         because there are derivatives in three directions. */
  vector<ColMajorMatrix<passivedouble> > derLagBasisIntLGL;  /*!< \brief The values of the derivatives of the Lagrangian
                                                                         basis functions in the integration points for the
                                                                         LGL point distribution. It is a vector, because
                                                                         there are derivatives in three directions. */

  /*!
   * \brief Function, which computes the values of the derivatives of the Lagrangian
   *        basis functions of a prism in the integration points for the given
   *        location of the DOFs.
   * \param[in]  rTriangleDOFs - Vector, which contains the parametric r-locations of the DOFs
   *                             of the base triangle.
   * \param[in]  sTriangleDOFs - Vector, which contains the parametric s-locations of the DOFs
   *                             of the base triangle.
   * \param[in]  rLineDOFs     - Vector, which contains the parametric r-locations of the DOFs
   *                             of the line normal to the triangles.
   * \param[out] derLag        - Matrix, which contains the values of derivatives of all the
   *                             Lagrangian basis functions in all the integration points.
   */
  void DerLagBasisIntPointsPrism(const vector<passivedouble>            &rTriangleDOFs,
                                 const vector<passivedouble>            &sTriangleDOFs,
                                 const vector<passivedouble>            &rLineDOFs,
                                 vector<ColMajorMatrix<passivedouble> > &derLag);

  /*!
   * \brief Function, which computes the values of the Lagrangian basis functions
   *        of a prism in the integration points for the given location of the DOFs.
   * \param[in]  rTriangleDOFs - Vector, which contains the parametric r-locations of the DOFs
   *                             of the base triangle.
   * \param[in]  sTriangleDOFs - Vector, which contains the parametric s-locations of the DOFs
   *                             of the base triangle.
   * \param[in]  rLineDOFs     - Vector, which contains the parametric r-locations of the DOFs
   *                             of the line normal to the triangles.
   * \param[out] lag           - Matrix, which contains the values of all the Lagrangian
   *                             basis functions in all the integration points.
   */
  void LagBasisIntPointsPrism(const vector<passivedouble>   &rTriangleDOFs,
                              const vector<passivedouble>   &sTriangleDOFs,
                              const vector<passivedouble>   &rLineDOFs,
                              ColMajorMatrix<passivedouble> &lag);

  /*!
   * \brief Function, which computes the location of all the DOFs of the prim
   *        from the given coordinates of the triangle and line.
   * \param[in]  rTriangleDOFs - Vector, which contains the parametric r-locations of the DOFs
   *                             of the base triangle.
   * \param[in]  sTriangleDOFs - Vector, which contains the parametric s-locations of the DOFs
   *                             of the base triangle.
   * \param[in]  rLineDOFs     - Vector, which contains the parametric r-locations of the DOFs
   *                             of the line normal to the triangles.
   * \param[out] rDOFs         - Parametric r-coordinates of the DOFs of the prism.
   * \param[out] sDOFs         - Parametric s-coordinates of the DOFs of the prism.
   * \param[out] tDOFs         - Parametric t-coordinates of the DOFs of the prism.
   */
  void LocationAllDOFsPrism(const vector<passivedouble>   &rTriangleDOFs,
                            const vector<passivedouble>   &sTriangleDOFs,
                            const vector<passivedouble>   &rLineDOFs,
                            vector<passivedouble>         &rDOFs,
                            vector<passivedouble>         &sDOFs,
                            vector<passivedouble>         &tDOFs);

  /*!
   * \brief Function, which computes the gradients of the Vandermonde matrix for a standard prism.
   * \param[in]  r   - Parametric r-coordinates for which the Vandermonde matrix must be computed.
   * \param[in]  s   - Parametric s-coordinates for which the Vandermonde matrix must be computed
   * \param[in]  t   - Parametric t-coordinates for which the Vandermonde matrix must be computed
   * \param[out] VDr - Matrix to store the derivative in r-direction of the Vandermonde matrix
   *                   in all r,s,t-locations.
   * \param[out] VDs - Matrix to store the derivative in s-direction of the Vandermonde matrix
   *                   in all r,s,t-locations.
   * \param[out] VDt - Matrix to store the derivative in t-direction of the Vandermonde matrix
   *                   in all r,s,t-locations.
   */
  void GradVandermondePrism(const vector<passivedouble>   &r,
                            const vector<passivedouble>   &s,
                            const vector<passivedouble>   &t,
                            ColMajorMatrix<passivedouble> &VDr,
                            ColMajorMatrix<passivedouble> &VDs,
                            ColMajorMatrix<passivedouble> &VDt);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard prism.
   * \param[in]  r  - Parametric r-coordinates for which the Vandermonde matrix must be computed.
   * \param[in]  s  - Parametric s-coordinates for which the Vandermonde matrix must be computed
   * \param[in]  t  - Parametric t-coordinates for which the Vandermonde matrix must be computed
   * \param[out] V  - Matrix to store the Vandermonde matrix in all r,s,t-locations.
   */
  void VandermondePrism(const vector<passivedouble>   &r,
                        const vector<passivedouble>   &s,
                        const vector<passivedouble>   &t,
                        ColMajorMatrix<passivedouble> &V);

  /*!
   * \brief Function, which creates the local grid connectivities of the faces
   *        of the volume element.
   */
  void LocalGridConnFaces(void);
};
