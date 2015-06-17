/*!
 * \file interpolation_structure.hpp
 * \brief Headers of the main subroutines used by SU2_FSI.
 *        The subroutines and functions are in the <i>interpolation_structure.cpp</i> file.
 * \author H. Kline
 * \version 3.2.9 "eagle"
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

#ifdef HAVE_MPI
  #include "mpi.h"
#endif

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
  unsigned short nZone;
  double ***Data; /*!\brief container for some data to be interpolated */
public:
  CGeometry*** Geometry; /*! \brief Vector which stores n zones of geometry. */

  /*!
   * \brief Constructor of the class.
   */
  CInterpolator(void);

  /*!
 * \brief Constructor of the class.
 */
  CInterpolator(CGeometry ***geometry_container, CConfig **config,  unsigned short iZone_0,unsigned short iZone_1, unsigned short nZone);

  /*!
   * \brief Destructor of the class.
   */
  ~CInterpolator(void);

  /*!
     * \brief initialize the Data structure to the appropriate size.
     */
  void InitializeData(unsigned short iZone_0, unsigned short iZone_1, unsigned short nVar);

  /*!
   * \brief interpolate Data from one mesh to another
   */
  void Interpolate_Data(unsigned short iZone_0, unsigned short iZone_1, CConfig **config);

  /*!
   * \brief interpolate deformations from one mesh to another
   */
  void Interpolate_Deformation(unsigned short iZone_0, unsigned short iZone_1, CConfig **config);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   */
  void Set_TransferCoeff(unsigned short iZone_0, unsigned short iZone_1, CConfig **config);


  /*!
   * \brief Return the value of the Data at the specified zone, point, and dimension.
   */
  double GetData(unsigned short iZone, unsigned long iPoint, unsigned short iDim);

  /*!
   * \brief Return the value of the Data vector at the specified zone and point.
   */
  double* GetData(unsigned short iZone, unsigned long iPoint);

  /*!
   * \brief Set the value of the Data at the specified zone, point, and dimension.
   */
  void SetData(unsigned short iZone, unsigned long iPoint, unsigned short iDim, double val);




};


class CNearestNeighbor : public CInterpolator {
public:

  /*!
   * \brief Constructor of the class.
   */
  CNearestNeighbor(CGeometry ***geometry_container, CConfig **config,  unsigned short iZone_0,unsigned short iZone_1,unsigned short nZone);

  /*!
   * \brief Destructor of the class.
   */
  ~CNearestNeighbor(void);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   */
  void Set_TransferCoeff(unsigned short iZone_0, unsigned short iZone_1, CConfig **config);

};
