/*!
 * \file interpolation_structure.hpp
 * \brief Headers of the main subroutines used by SU2_FSI.
 *        The subroutines and functions are in the <i>interpolation_structure.cpp</i> file.
 * \author H. Kline
 * \version 4.0.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../../Common/include/mpi_structure.hpp"

#include <iostream>
#include <cmath>
#include <cstdlib>

#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "vector_structure.hpp"

using namespace std;




/*!
 * \class CInterpolator
 * \brief Main class for defining the interpolator, it requires
 * a child class for each particular interpolation method
 * \author H. Kline
 * \version 3.2.9 "eagle"
 */
class CInterpolator {
protected:
  unsigned int nZone;
  unsigned short nVar;
  su2double ***Data; /*!\brief container for some data to be interpolated */
public:
  CGeometry*** Geometry; /*! \brief Vector which stores n zones of geometry. */

  /*!
   * \brief Constructor of the class.
   */
  CInterpolator(void);

  /*!
 * \brief Constructor of the class.
 */
  CInterpolator(CGeometry ***geometry_container, CConfig **config,  unsigned int* Zones, unsigned int nZone);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CInterpolator(void);

  /*!
     * \brief initialize the Data structure to the appropriate size.
     */
  void InitializeData(unsigned int* Zones, unsigned short val_nVar);

  /*!
   * \brief interpolate Data from one mesh to another.
   * The data for zone 0 will be overwritten. transfer coefficients must be defined with Set_TransferCoeff.
   * \param[in] iZone_0 - zone to recieve interpolated data
   * \param[in] config
   */
  void Interpolate_Data(unsigned int iZone,  CConfig **config);

  /*!
   * \brief interpolate deformations from one mesh to another.
   * Uses information stored by the geometry class, updates values in VarCoord of iZone_0. Set_TransferCoeff must be run first.
   * \param[in] iZone_0 - zone to recieve interpolated data.
   * \param[in] config
   */
  void Interpolate_Deformation(unsigned int iZone, CConfig **config);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   * \param[in] Zones - list of zones to set up interpolation for. This method must be overwritten in the child classes.
   * \param[in] config
   */
  virtual void Set_TransferCoeff(unsigned int* Zones, CConfig **config)=0;


  /*!
   * \brief Return the value of the Data at the specified zone, point, and dimension.
   * \param[in] iZone - zone index
   * \param[in] iPoint - point index
   * \param[in[ iDim - index of the data
   */
  su2double GetData(unsigned int iZone, unsigned long iPoint, unsigned short iVar);

  /*!
   * \brief Return the pointer to the Data vector at the specified zone and point.
   */
  su2double* GetData(unsigned int iZone, unsigned long iPoint);

  /*!
   * \brief Set the value of the Data at the specified zone, point, and index.
   */
  void SetData(unsigned int iZone, unsigned long iPoint, unsigned short iVar, su2double val);




};

/*!
 * \brief Nearest Neighbor interpolation
 */
class CNearestNeighbor : public CInterpolator {
public:

  /*!
   * \brief Constructor of the class.
   */
  CNearestNeighbor(CGeometry ***geometry_container, CConfig **config,  unsigned int* Zones,unsigned int nZone);

  /*!
   * \brief Destructor of the class.
   */
  ~CNearestNeighbor(void);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   */
  void Set_TransferCoeff(unsigned int* Zones, CConfig **config);

};

/*!
 * \brief Consistent and Conservative interpolation
 */
class CConsistConserve : public CInterpolator {
public:

  /*!
   * \brief Constructor of the class.
   */
  CConsistConserve(CGeometry ***geometry_container, CConfig **config,  unsigned int* Zones,unsigned int nZone);

  /*!
   * \brief Destructor of the class.
   */
  ~CConsistConserve(void);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   * \param[in] Zones - list of zones to use for interpolation
   * \param[in] config -
   */
  void Set_TransferCoeff(unsigned int* Zones, CConfig **config);

  /*!
   * \brief Calculate the isoparametric representation of point iVertex in marker iZone_0 by nodes of element donor_elem in marker jMarker of zone iZone_1.
   * \param[out] isoparams - isoparametric coefficients. Must be allocated to size nNodes ahead of time.
   * \param[in] iZone_0 - zone index of the point being interpolated
   * \param[in] iMarker - marker index of the point being interpolated
   * \param[in] iVertex - vertex index of the point being interpolated.
   * \param[in] nDim - the dimension of the coordinates.
   * \param[in] iZone_1 - zone index of the element to use for interpolation
   * \param[in] jMarker - marker index of the element to use for interpolation
   * \param[in] donor_elem - element index of the element to use for interpolation
   * \param[in[ nDonorPoints - number of donor points in the element.
   */
  void Isoparametric(su2double* isoparams, unsigned int iZone_0, unsigned short iMarker, unsigned long iVertex, unsigned int nDim, unsigned int iZone_1, unsigned short jMarker, long donor_elem,  unsigned short iFace, unsigned int nDonorPoints);

};
