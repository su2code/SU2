/*!
 * \file vector_structure.inl
 * \brief inline subroutines of the <i>vector_structure.hpp</i> file.
 * \author F. Palacios, J. Hicken
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
inline void CSysVector<ScalarType>::SetValZero(void) { 
  for (unsigned long i = 0; i < nElm; i++)
		vec_val[i] = 0.0;
}

template<class ScalarType>
inline unsigned long CSysVector<ScalarType>::GetLocSize() const { return nElm; }

template<class ScalarType>
unsigned long CSysVector<ScalarType>::GetNElmDomain() const { return nElmDomain; }

template<class ScalarType>
inline unsigned long CSysVector<ScalarType>::GetSize() const {
#ifdef HAVE_MPI
  return nElmGlobal;
#else
  return (unsigned long)nElm;
#endif
}

template<class ScalarType>
inline unsigned short CSysVector<ScalarType>::GetNVar() const { return nVar; }

template<class ScalarType>
inline unsigned long CSysVector<ScalarType>::GetNBlk() const { return nBlk; }

template<class ScalarType>
inline unsigned long CSysVector<ScalarType>::GetNBlkDomain() const { return nBlkDomain; }

template<class ScalarType>
inline ScalarType & CSysVector<ScalarType>::operator[](const unsigned long & i) { return vec_val[i]; }

template<class ScalarType>
inline const ScalarType & CSysVector<ScalarType>::operator[](const unsigned long & i) const { return vec_val[i]; }

template<class ScalarType>
template<class T>
void CSysVector<ScalarType>::PassiveCopy(const CSysVector<T>& other) {

  /*--- This is a method and not the overload of an operator to make sure who
   calls it knows the consequence to the derivative information (lost) ---*/

  /*--- check if self-assignment, otherwise perform deep copy ---*/
  if ((const void*)this == (const void*)&other) return;

  /*--- determine if (re-)allocation is needed ---*/
  if (nElm != other.GetLocSize() && vec_val != NULL) {
    delete [] vec_val;
    vec_val = NULL;
  }

  /*--- copy ---*/
  nElm = other.GetLocSize();
  nElmDomain = other.GetNElmDomain();
  nBlk = other.GetNBlk();
	nBlkDomain = other.GetNBlkDomain();
  nVar = other.GetNVar();

  if (vec_val == NULL)
    vec_val = new ScalarType[nElm];

  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = SU2_TYPE::GetValue(other[i]);

#ifdef HAVE_MPI
  nElmGlobal = other.GetSize();
#endif
}
