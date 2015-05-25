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
#include <ctime>

#include "geometry_structure.hpp"
#include "config_structure.hpp"

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
public:
  CGeometry** Geometry; /*! \brief Vector which stores n zones of geometry. */
  CSysTransferMatrix* TransferMatrix; /*! \brief Sparse matrix structure defining transfer from one mesh to another. */

  /*!
   * \brief Constructor of the class.
   */
  CInterpolator(void);

  /*!
 * \brief Constructor of the class.
 */
  CInterpolator(CGeometry **geometry_container, CConfig **config, unsigned short nZone);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CInterpolator(void);

  /*!
   * \brief interpolate forces from one mesh to another
   */
  void Interpolate_Force(unsigned short iZone_0, unsigned short iZone_1);

  /*!
   * \brief interpolate deformations from one mesh to another
   */
  void Interpolate_Deformation(unsigned short iZone_0, unsigned short iZone_1);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   */
  void Set_TransferMatrix(unsigned short iZone_0, unsigned short iZone_1);


};


class CNearestNeighbor : public CInterpolator {
public:
  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   */
  void Set_TransferMatrix(unsigned short iZone_0, unsigned short iZone_1);

};
