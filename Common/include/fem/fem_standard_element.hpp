/*!
 * \file fem_standard_element.hpp
 * \brief Headers of the main functions for the FEM standard elements.
 *        The functions are in the <i>fem_standard_element.cpp</i> file.
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

#include <iostream>
#include <vector>
#include <cstdlib>

#include "../CConfig.hpp"

using namespace std;

/*!
 * \class CFEMStandardElementBase
 * \brief Base class for a FEM standard element.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
class CFEMStandardElementBase {
protected:
  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
  unsigned short orderExact;   /*!< \brief Polynomial order that must be integrated exactly by the integration rule. */
  unsigned short nIntegration; /*!< \brief Number of points used in the numerical integration. */

public:
  /*!
  * \brief Constructor. Nothing to be done.
  */
  CFEMStandardElementBase() = default;

  /*!
  * \brief Destructor. Nothing to be done, because the vectors are deleted automatically.
  */
  virtual ~CFEMStandardElementBase() = default;

public:
  /*!
  * \brief Function, which makes available the type of the element.
  * \return  The type of the element using the VTK convention.
  */
  inline unsigned short GetVTK_Type(void) const {return VTK_Type;}

  /*!
  * \brief Static function, which makes available the number of DOFs for an element
           corresponding to the arguments.
  * \param[in] VTK_Type   - Type of the element using the VTK convention.
  * \param[in] nPoly      - Polynomial degree of the element.
  * \return The number of DOFs
  */
  static unsigned short GetNDOFsStatic(unsigned short VTK_Type,
                                       unsigned short nPoly);

protected:
  /*!
  * \brief Function, which computes the local connectivity of linear subelements of
           a line, which can be used for plotting.
  * \param[in]  nPoly    - Polynomial degree of the line.
  * \param[out] subConn  - The local subconnectivity of a line element.
  */
  void SubConnForPlottingLine(const unsigned short   nPoly,
                              vector<unsigned short> &subConn);

  /*!
  * \brief Function, which computes the local connectivity of linear subelements of
           a quadrilateral, which can be used for plotting.
  * \param[in]  nPoly    - Polynomial degree of the quadrilateral.
  * \param[out] subConn  - The local subconnectivity of a triangle element.
  */
  void SubConnForPlottingQuadrilateral(const unsigned short   nPoly,
                                       vector<unsigned short> &subConn);

  /*!
  * \brief Function, which computes the local connectivity of linear subelements of
           a triangle, which can be used for plotting.
  * \param[in]  nPoly    - Polynomial degree of the triangle.
  * \param[out] subConn  - The local subconnectivity of a triangle element.
  */
  void SubConnForPlottingTriangle(const unsigned short   nPoly,
                                  vector<unsigned short> &subConn);
};

/*!
 * \class CFEMStandardElement
 * \brief Class to define a FEM standard element.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
class CFEMStandardElement : public CFEMStandardElementBase {
private:

  unsigned short nPoly;        /*!< \brief Polynomial degree of the element. */
  unsigned short nDOFs;        /*!< \brief Number of DOFs of the element. */

  unsigned short VTK_Type1;    /*!< \brief VTK type for elements of type 1 in subConn1ForPlotting. */
  unsigned short VTK_Type2;    /*!< \brief VTK type for elements of type 2 in subConn2ForPlotting. */

public:
  /*!
   * \brief Function, which makes available the type of the element in subConn1ForPlotting.
   * \return  The type of the elements in subConn1ForPlotting using the VTK convention.
   */
  inline unsigned short GetVTK_Type1(void) const {return VTK_Type1;}

  /*!
   * \brief Function, which makes available the number of sub-elements of type 1 for plotting.
   * \return  The number of sub-elements of type 1 for plotting.
   */
  inline unsigned short GetNSubElemsType1(void) const {return 0;}

  /*!
   * \brief Function, which makes available the the connectivity of the linear elements of type 1 as a const pointer.
   * \return  The pointer to the local connectivity of the linear elements of type 1.
   */
  inline const unsigned short *GetSubConnType1(void) const {return NULL;}

  /*!
   * \brief Function, which makes available the type of the element in subConn2ForPlotting.
   * \return  The type of the elements in subConn2ForPlotting using the VTK convention.
   */
  inline unsigned short GetVTK_Type2(void) const {return VTK_Type2;}

  /*!
   * \brief Function, which makes available the number of sub-elements of type 2 for plotting.
   * \return  The number of sub-elements of type 2 for plotting.
   */
  inline unsigned short GetNSubElemsType2(void) const {return 0;}

  /*!
   * \brief Function, which makes available the the connectivity of the linear elements of type 2 as a const pointer.
   * \return  The pointer to the local connectivity of the linear elements of type 2.
   */
  inline const unsigned short *GetSubConnType2(void) const {return NULL;}

  /*!
   * \brief Function, which makes available the number of DOFs of a linear element, used for plotting.
   * \return  The number of DOFs of the linear elements.
   */
  unsigned short GetNDOFsPerSubElem(unsigned short val_VTK_Type) const;
};

/*!
 * \class CFEMStandardBoundaryFace
 * \brief Class to define a FEM standard boundary face.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
class CFEMStandardBoundaryFace : public CFEMStandardElementBase {

public:
  /*!
  * \brief Function, which makes available the number of DOFs of a linear subface, used
           for plotting, among others, plotting.
  * \return  The number of DOFs of a linear subfaces of the face.
  */
  inline unsigned short GetNDOFsPerSubFace(void) const {return 0;}

  /*!
  * \brief Function, which makes available the number of linear subfaces used
           for plotting, among others.
  * \return  The number of linear subfaces of the face.
  */
  inline unsigned short GetNSubFaces(void) const {return 0;}

  /*!
  * \brief Function, which makes available the the connectivity of the linear subfaces
           as a const pointer.
  * \return  The pointer to the local connectivity of the linear subfaces.
  */
  inline const unsigned short* GetSubFaceConn(void) const {return NULL;}
};
