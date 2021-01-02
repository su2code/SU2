/*!
 * \file CGemmStandard.cpp
 * \brief Functions for the class CGemmStandard.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
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

#include "../../include/fem/CGemmStandard.hpp"
#include "../../include/fem/CFEMStandardElementBase.hpp"

/*----------------------------------------------------------------------------------*/
/*                  Public member functions of CGemmStandard.                       */
/*----------------------------------------------------------------------------------*/

CGemmStandard::CGemmStandard(const int val_M, const int val_N, const int val_K)
  : CGemmBase() {

  /*--- Copy the arguments into the member variables. ---*/
  M = val_M;
  N = val_N;
  K = val_K;

  /*--- Create the padded value of M. ---*/
  MP = CFEMStandardElementBase::PaddedValue(M);
}

CGemmStandard::~CGemmStandard(){}
