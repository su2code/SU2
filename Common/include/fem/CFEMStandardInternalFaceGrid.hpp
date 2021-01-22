/*!
 * \file CFEMStandardInternalFaceGrid.hpp
 * \brief Class for the FEM standard element for internal matching faces
 *        for the grid.
 *        The functions are in the <i>CFEMStandardInternalFaceGrid.cpp</i> file.
 * \author E. van der Weide
 * \version 7.1.0 "Blackbird"
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

#include "CFEMStandardElementBase.hpp"

/*!
 * \class CFEMStandardInternalFaceGrid
 * \brief Class which defines the variables and methods for the standard
 *        standard element for an internal matching face for the grid.
 * \author E. van der Weide
 * \version 7.1.0 "Blackbird"
 */
class CFEMStandardInternalFaceGrid final {
private:
  CFEMStandardElementBase *elem0;  /*!< \brief Standard element on side 0 of the internal matching face. */
  CFEMStandardElementBase *elem1;  /*!< \brief Standard element on side 1 of the internal matching face. */

public:
  /*-----------------------------------------------------------------------------------*/
  /*---                     Constructors and destructors.                           ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardInternalFaceGrid() = delete;

  /*!
   * \overload
   * \param[in] val_elem0 - Standard surface elements for side 0.
   * \param[in] val_elem1 - Standard surface elements for side 1
   */
  CFEMStandardInternalFaceGrid(CFEMStandardElementBase *val_elem0,
                               CFEMStandardElementBase *val_elem1);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  ~CFEMStandardInternalFaceGrid() = default;

  /*-----------------------------------------------------------------------------------*/
  /*---                  Inline public member functions.                            ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which makes available the number of total integration points of the element.
   * \return  The number of total integration points.
   */
  inline unsigned short GetNIntegration(void) const {return elem0->GetNIntegration();}

  /*!
   * \brief Function, which makes available the padded number of total integration points of the element.
   * \return  The padded number of total integration points.
   */
  inline unsigned short GetNIntegrationPad(void) const {return elem0->GetNIntegrationPad();}

  /*-----------------------------------------------------------------------------------*/
  /*---                     Public member functions.                                ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which determines the coordinates of the integration points of the
   *        face from the coordinates of the element on side 0. This is OK, because the
   *        the polynomial order of the grid of the element on side 1 is not higher than
   *        the polynomial order of the grid of the element  on side 0, 
   * \param[in]  coorGridDOFsVol   - The coordinates of the grid DOFs of the element on side 0.
   * \param[out] coorIntPointsFace - The coordinates of the integration points on the face.
   */
  void CoorIntPoints(ColMajorMatrix<su2double> &coorGridDOFsVol,
                     ColMajorMatrix<su2double> &coorIntPointsFace);

  /*!
   * \brief Function, which determines the coordinates of the integration points of the
   *        face from the coordinates of the element on side 1. This function is only
   *        used for debugging purposes.
   * \param[in]  coorGridDOFsVol   - The coordinates of the grid DOFs of the element on side 0.
   * \param[out] coorIntPointsFace - The coordinates of the integration points on the face.
   */
  void CoorIntPointsFromSide1(ColMajorMatrix<su2double> &coorGridDOFsVol,
                              ColMajorMatrix<su2double> &coorIntPointsFace);

  /*!
   * \brief Function, which computes the metric terms in the surface integration points.
   * \param[in]  matCoorElem0     - Matrix that contains the coordinates of the grid DOFs
   *                                of the element on side 0.
   * \param[in]  matCoorElem1     - Matrix that contains the coordinates of the grid DOFs
   *                                of the element on side 1.
   * \param[out] JacobiansFace    - Jacobians of the integration points of the face.
   * \param[out] normalsFace      - Unit normals in the integration points of the face.
   * \param[out] metricTermsSide0 - The metric terms drdx, drdy, etc. in the integration
   *                                points on side 0 of the face.
   * \param[out] metricTermsSide1 - The metric terms drdx, drdy, etc. in the integration
   *                                points on side 1 of the face.
   */
   void MetricTermsSurfaceIntPoints(ColMajorMatrix<su2double>          &matCoorElem0,
                                    ColMajorMatrix<su2double>          &matCoorElem1,
                                    su2activevector                    &JacobiansFace,
                                    ColMajorMatrix<su2double>          &normalsFace,
                                    vector<ColMajorMatrix<su2double> > &metricTermsSide0,
                                    vector<ColMajorMatrix<su2double> > &metricTermsSide1);
private:

};
