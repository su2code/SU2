/*!
 * \file adt_structure.inl
 * \brief In-Line subroutines of the <i>adt_structure.hpp</i> file.
 * \author E. van der Weide
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

inline CADTComparePointClass::~CADTComparePointClass() {}

inline bool CADTComparePointClass::operator()(const unsigned long p0,
                                              const unsigned long p1) const {
  return pointCoor[nDim*p0+splitDirection] < pointCoor[nDim*p1+splitDirection];
}

inline CBBoxTargetClass::CBBoxTargetClass() {}

inline CBBoxTargetClass::~CBBoxTargetClass() {}

inline CBBoxTargetClass::CBBoxTargetClass(const CBBoxTargetClass &other) {Copy(other);}

inline CBBoxTargetClass& CBBoxTargetClass::operator=(const CBBoxTargetClass &other) {Copy(other); return (*this);}

inline CADTNodeClass::CADTNodeClass() {}

inline CADTNodeClass::~CADTNodeClass() {}

inline CADTNodeClass::CADTNodeClass(const CADTNodeClass &other) {Copy(other);}

inline CADTNodeClass& CADTNodeClass::operator=(const CADTNodeClass &other) {Copy(other); return (*this);}

inline CADTBaseClass::CADTBaseClass() {}

inline CADTBaseClass::~CADTBaseClass() {}

inline bool CADTBaseClass::IsEmpty(void) const { return isEmpty;}

inline CADTPointsOnlyClass::~CADTPointsOnlyClass() {}

inline CADTElemClass::~CADTElemClass() {}

