/*!
 * \file CFEMStandardVolumeHexGrid.hpp
 * \brief Class for the FEM volume hexahedron standard element for the grid.
 *        The functions are in the <i>CFEMStandardVolumeHexGrid.cpp</i> file.
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

#include "CFEMStandardHex.hpp"

/*!
 * \class CFEMStandardVolumeHexGrid
 * \brief Class which defines the variables and methods for the volume
 *        hexahedron standard element for the grid.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
class CFEMStandardVolumeHexGrid final: public CFEMStandardHex {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardVolumeHexGrid() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardVolumeHexGrid(const unsigned short val_nPoly,
                            const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  virtual ~CFEMStandardVolumeHexGrid() = default;

  /*!
   * \brief Function to determine whether the point distribution of the element is LGL.
   * \param[in] matCoor - Matrix that contains the coordinates of the grid DOFs.
   * \param[in] ldb     - Leading dimension of matCoor (gemm convention).
   * \return True if the coordinates are LGL and false otherwise.
   */
  bool CoordinatesAreLGL(const ColMajorMatrix<su2double> &matCoor,
                         const unsigned short            ldb) const override;

  /*!
   * \brief Function, which computes the data and/or derivatives in the
   *        integration points from the known data in the DOFs.
   * \param[in]  matB    - Matrix that contains the input data.
   * \param[in]  ldb     - Leading dimension of matB (gemm convention).
   * \param[in]  ldc     - Leading dimension of matC (gemm convention).
   * \param[in]  n       - Second dimension of matB and matC (gemm convention).
   * \param[out] matC    - Result of the multiplication C = A*B.
   * \param[out] matDerC - Result of the multiplication CDer = ADer*B.
   * \param[in]  config  - Pointer to the configuration. Used for the timings.
   */
  void DataIntegrationPoints(const ColMajorMatrix<su2double>    &matB,
                             const unsigned short               ldb,
                             const unsigned short               ldc,
                             const unsigned short               n,
                             ColMajorMatrix<su2double>          *matC,
                             vector<ColMajorMatrix<su2double> > *matDerC,
                             const CConfig                      *config) const override;

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
  unsigned short GetVTK_TypeFace(unsigned short ind) const override {return QUADRILATERAL;}

private:

};
