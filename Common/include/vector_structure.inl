/*!
 * \file vector_structure.inl
 * \brief inline subroutines of the <i>vector_structure.hpp</i> file.
 * \author F. Palacios, J. Hicken
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

template<class CalcType>
inline void TCSysVector<CalcType>::SetValZero(void) { 
  for (unsigned long i = 0; i < nElm; i++)
		vec_val[i] = 0.0;
}

template<class CalcType>
inline unsigned long TCSysVector<CalcType>::GetLocSize() const { return nElm; }

template<class CalcType>
inline unsigned long TCSysVector<CalcType>::GetSize() const {
#ifdef HAVE_MPI
  return nElmGlobal;
#else
  return (unsigned long)nElm;
#endif
}

template<class CalcType>
inline unsigned short TCSysVector<CalcType>::GetNVar() const { return nVar; }

template<class CalcType>
inline unsigned long TCSysVector<CalcType>::GetNBlk() const { return nBlk; }

template<class CalcType>
inline unsigned long TCSysVector<CalcType>::GetNBlkDomain() const { return nBlkDomain; }

template<class CalcType>
inline CalcType & TCSysVector<CalcType>::operator[](const unsigned long & i) { return vec_val[i]; }

template<class CalcType>
inline const CalcType & TCSysVector<CalcType>::operator[](const unsigned long & i) const { return vec_val[i]; }
