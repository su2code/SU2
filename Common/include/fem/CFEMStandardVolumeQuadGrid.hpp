/*!
 * \file CFEMStandardVolumeQuadGrid.hpp
 * \brief Class for the FEM volume quadrilateral standard element for the grid.
 *        The functions are in the <i>CFEMStandardVolumeQuadGrid.cpp</i> file.
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

#include "CFEMStandardQuad.hpp"

/*!
 * \class CFEMStandardVolumeQuadGrid
 * \brief Class which defines the variables and methods for the volume
 *        quadrilateral standard element for the grid.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
class CFEMStandardVolumeQuadGrid final: public CFEMStandardQuad {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardVolumeQuadGrid() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardVolumeQuadGrid(const unsigned short val_nPoly,
                             const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  ~CFEMStandardVolumeQuadGrid() = default;

  /*!
   * \brief Function, which computes the derivatives of the coordinates in the volume
   *        integration points.
   * \param[in]  LGLDistribution - Whether or not the LGL node distribution must be used.
   * \param[in]  matCoor         - Matrix that contains the coordinates of the grid DOFs.
   * \param[out] matDerCoor      - Vector of matrices to store the derivatives of the coordinates.
   */
  void DerivativesCoorVolumeIntPoints(const bool                         LGLDistribution,
                                      const ColMajorMatrix<su2double>    &matCoor,
                                      vector<ColMajorMatrix<su2double> > &matDerCoor) const override;

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
};
