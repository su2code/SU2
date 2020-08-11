/*!
 * \file CFEMStandardVolumeTriGrid.cpp
 * \brief Functions for the class CFEMStandardVolumeTriGrid.
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

#include "../../include/fem/CFEMStandardVolumeTriGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*          Public member functions of CFEMStandardVolumeTriGrid.                   */
/*----------------------------------------------------------------------------------*/

CFEMStandardVolumeTriGrid::CFEMStandardVolumeTriGrid(const unsigned short val_nPoly,
                                                     const unsigned short val_orderExact) {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = TRIANGLE;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;
}
