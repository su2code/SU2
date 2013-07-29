/*!
 * \file vector_structure.inl
 * \brief inline subroutines of the <i>vector_structure.hpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
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

inline double & CSysVector::operator[](const unsigned int & i) { return vec_val[i]; }

inline const double & CSysVector::operator[](const unsigned int & i) const { return vec_val[i]; }
