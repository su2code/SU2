/*!
 * \file fem_geometry_structure.inl
 * \brief In-Line subroutines of the <i>fem_geometry_structure.hpp</i> file.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
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
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

inline CPointCompare::CPointCompare(void) { }

inline CPointCompare::~CPointCompare(void) { }

inline CPointCompare::CPointCompare(const CPointCompare &other) { Copy(other); }

inline CPointCompare& CPointCompare::operator=(const CPointCompare &other) { Copy(other); return (*this); }

inline CVolumeElementFEM::CVolumeElementFEM(void) { }

inline CVolumeElementFEM::~CVolumeElementFEM(void) { }

inline CPointFEM::CPointFEM(void) { coor[0] = coor[1] = coor[2] = 0.0; }

inline CPointFEM::~CPointFEM(void) { }

inline CPointFEM::CPointFEM(const CPointFEM &other) { Copy(other); }

inline CPointFEM& CPointFEM::operator=(const CPointFEM &other) { Copy(other); return (*this); }

inline CSurfaceElementFEM::CSurfaceElementFEM(void) { indStandardElement = -1; }

inline CSurfaceElementFEM::~CSurfaceElementFEM(void) { }

inline CSurfaceElementFEM::CSurfaceElementFEM(const CSurfaceElementFEM &other) { Copy(other); }

inline CSurfaceElementFEM& CSurfaceElementFEM::operator=(const CSurfaceElementFEM &other) { Copy(other); return (*this); }

inline bool CSurfaceElementFEM::operator< (const CSurfaceElementFEM &other) const { return boundElemIDGlobal < other.boundElemIDGlobal; }

inline CBoundaryFEM::CBoundaryFEM(void) { }

inline CBoundaryFEM::~CBoundaryFEM(void) { }

inline CMeshFEM::CMeshFEM(void) { }

inline CMeshFEM::~CMeshFEM(void) { }

inline CMeshFEM_DG::CMeshFEM_DG(void) { }

inline CMeshFEM_DG::~CMeshFEM_DG(void) { }
