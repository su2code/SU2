/*!
 * \file linear_solvers_structure.inl
 * \brief inline subroutines of the <i>linear_solvers_structure.hpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.3
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

inline unsigned int CSysVector::GetLocSize() const { return nElm; }

inline unsigned long CSysVector::GetSize() const {
#ifndef NO_MPI
  return nElmGlobal;
#else
  return (unsigned long)nElm;
#endif
}

inline unsigned short CSysVector::GetNVar() const { return nVar; }

inline unsigned int CSysVector::GetNBlk() const { return nBlk; }

inline unsigned int CSysVector::GetNBlkDomain() const { return nBlkDomain; }

inline double & CSysVector::operator[](const unsigned int & i) {
#ifdef DEBUG
  /*--- check for invalid index ---*/
  if ( (i < 0) || (i >= nElm) ) {
    cerr << "CSysVector::operator[](unsigned int): " 
	 << "invalid index value: i = " << i << endl;
    throw(-1);
  }
#endif
  return vec_val[i];
}

inline const double & CSysVector::operator[](const unsigned int & i) const {
#ifdef DEBUG
  /*--- check for invalid index ---*/
  if ( (i < 0) || (i >= nElm) ) {
    cerr << "CSysVector::operator[](unsigned int) const: " 
	 << "invalid index value: i = " << i << endl;
    throw(-1);
  }
#endif
  return vec_val[i];
}

inline double CSysSolve::sign(const double & x, const double & y) const {
  if (y == 0.0)
    return 0.0;
  else {
    return (y < 0 ? -fabs(x) : fabs(x));
  }
}
