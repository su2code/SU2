/*!
 * \file linear_solvers_structure.inl
 * \brief inline subroutines of the <i>linear_solvers_structure.hpp</i> file.
 * \author J. Hicken, F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

template<class ScalarType>
inline ScalarType CSysSolve<ScalarType>::Sign(const ScalarType & x, const ScalarType & y) const {
  if (y == 0.0)
    return 0.0;
  else {
//    return (y < 0 ? -fabs(x) : fabs(x));
    if (y < 0) return -fabs(x);
    else return fabs(x);
  }
}

template<class ScalarType>
inline ScalarType CSysSolve<ScalarType>::GetResidual(void) const { return Residual; }

template<class ScalarType>
void CSysSolve<ScalarType>::HandleTemporariesIn(CSysVector<su2double> & LinSysRes, CSysVector<su2double> & LinSysSol) {}

template<class ScalarType>
void CSysSolve<ScalarType>::HandleTemporariesOut(CSysVector<su2double> & LinSysSol) {}
