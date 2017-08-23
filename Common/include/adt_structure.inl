/*!
 * \file adt_structure.inl
 * \brief In-Line subroutines of the <i>adt_structure.hpp</i> file.
 * \author E. van der Weide
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
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

inline su2_adtComparePointClass::~su2_adtComparePointClass() {}

inline bool su2_adtComparePointClass::operator()(const unsigned long p0,
                                                 const unsigned long p1) const {
  return pointCoor[nDim*p0+splitDirection] < pointCoor[nDim*p1+splitDirection];
}

inline su2_adtNodeClass::su2_adtNodeClass() {}

inline su2_adtNodeClass::~su2_adtNodeClass() {}

inline su2_adtNodeClass::su2_adtNodeClass(const su2_adtNodeClass &other) {Copy(other);}

inline su2_adtNodeClass& su2_adtNodeClass::operator=(const su2_adtNodeClass &other) {Copy(other); return (*this);}

inline su2_adtBaseClass::su2_adtBaseClass() {}

inline su2_adtBaseClass::~su2_adtBaseClass() {}

inline bool su2_adtBaseClass::IsEmpty(void) const { return isEmpty;}

inline su2_adtPointsOnlyClass::~su2_adtPointsOnlyClass() {}
