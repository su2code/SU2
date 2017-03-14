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

inline long3T::long3T(){long0 = long1 = long2 = 0;}

inline long3T::long3T(const long a, const long b, const long c){long0 = a; long1 = b; long2 = c;}

inline long3T::~long3T(){}

inline long3T::long3T(const long3T &other){Copy(other);}

inline long3T& long3T::operator=(const long3T &other){Copy(other); return (*this);}

inline CReorderElementClass::CReorderElementClass(const unsigned long  val_GlobalElemID,
                                                  const unsigned short val_TimeLevel,
                                                  const bool           val_CommSolution)
 {globalElemID = val_GlobalElemID; timeLevel = val_TimeLevel; commSolution = val_CommSolution;}

inline CReorderElementClass::~CReorderElementClass(void) { }

inline CReorderElementClass::CReorderElementClass(const CReorderElementClass &other) { Copy(other); }

inline CReorderElementClass& CReorderElementClass::operator=(const CReorderElementClass &other) { Copy(other); return (*this); }

inline bool CReorderElementClass::GetCommSolution(void) { return commSolution; }

inline unsigned long CReorderElementClass::GetGlobalElemID(void) { return globalElemID; }

inline unsigned short CReorderElementClass::GetTimeLevel(void) { return timeLevel; }

inline SortFacesClass::SortFacesClass(unsigned long val_nVolElemOwned,
                                      unsigned long val_nVolElemTot)
 {nVolElemOwned = val_nVolElemOwned; nVolElemTot = val_nVolElemTot;}

inline SortFacesClass::~SortFacesClass(void) { }

inline CVolumeElementFEM::CVolumeElementFEM(void) { }

inline CVolumeElementFEM::~CVolumeElementFEM(void) { }

inline CPointFEM::CPointFEM(void) { coor[0] = coor[1] = coor[2] = 0.0; }

inline CPointFEM::~CPointFEM(void) { }

inline CPointFEM::CPointFEM(const CPointFEM &other) { Copy(other); }

inline CPointFEM& CPointFEM::operator=(const CPointFEM &other) { Copy(other); return (*this); }

inline CInternalFaceElementFEM::CInternalFaceElementFEM(void) { }

inline CInternalFaceElementFEM::~CInternalFaceElementFEM(void) { }

inline CSurfaceElementFEM::CSurfaceElementFEM(void) {indStandardElement = -1;}

inline CSurfaceElementFEM::~CSurfaceElementFEM(void) { }

inline CSurfaceElementFEM::CSurfaceElementFEM(const CSurfaceElementFEM &other) { Copy(other); }

inline CSurfaceElementFEM& CSurfaceElementFEM::operator=(const CSurfaceElementFEM &other) { Copy(other); return (*this); }

inline bool CSurfaceElementFEM::operator< (const CSurfaceElementFEM &other) const { return boundElemIDGlobal < other.boundElemIDGlobal; }

inline CBoundaryFEM::CBoundaryFEM(void) { }

inline CBoundaryFEM::~CBoundaryFEM(void) { }

inline CMeshFEM::CMeshFEM(void) { }

inline CMeshFEM::~CMeshFEM(void) { }

inline CBoundaryFEM* CMeshFEM::GetBoundaries(void) {return boundaries.data();}

inline CPointFEM* CMeshFEM::GetMeshPoints(void) {return meshPoints.data();}

inline unsigned long CMeshFEM::GetNMeshPoints(void) {return meshPoints.size();}

inline unsigned long CMeshFEM::GetNVolElemOwned(void) {return nVolElemOwned;}

inline unsigned long CMeshFEM::GetNVolElemTot(void) {return nVolElemTot;}

inline CVolumeElementFEM* CMeshFEM::GetVolElem(void) {return volElem.data();}

inline unsigned short CMeshFEM::GetNStandardBoundaryFacesSol(void) {return standardBoundaryFacesSol.size();}

inline FEMStandardBoundaryFaceClass* CMeshFEM::GetStandardBoundaryFacesSol(void) {return standardBoundaryFacesSol.data();}

inline const vector<int>& CMeshFEM::GetRanksRecv(void) const {return ranksRecv;}

inline const vector<int>& CMeshFEM::GetRanksSend(void) const {return ranksSend;}

inline const vector<vector<unsigned long> >& CMeshFEM::GetEntitiesRecv(void) const {return entitiesRecv;}

inline const vector<vector<unsigned long> >& CMeshFEM::GetEntitiesSend(void) const {return entitiesSend;}

inline vector<unsigned short> CMeshFEM::GetRotPerMarkers(void) const {return rotPerMarkers;}

inline vector<vector<unsigned long> > CMeshFEM::GetRotPerHalos(void) const {return rotPerHalos;}

inline CMeshFEM_DG::CMeshFEM_DG(void) { }

inline CMeshFEM_DG::~CMeshFEM_DG(void) { }

inline su2double* CMeshFEM_DG::GetLagrangianBeginTimeIntervalADER_DG(void) {return LagrangianBeginTimeIntervalADER_DG.data();}

inline su2double* CMeshFEM_DG::GetTimeInterpolDOFToIntegrationADER_DG(void) {return timeInterpolDOFToIntegrationADER_DG.data();}

inline unsigned long CMeshFEM_DG::GetNMatchingFacesWithHaloElem(void) {return nMatchingFacesWithHaloElem;}

inline unsigned long CMeshFEM_DG::GetNMatchingFaces(void) {return matchingFaces.size();}

inline CInternalFaceElementFEM* CMeshFEM_DG::GetMatchingFaces(void) {return matchingFaces.data();}

inline unsigned short CMeshFEM_DG::GetNStandardElementsSol(void) {return standardElementsSol.size();}

inline FEMStandardElementClass* CMeshFEM_DG::GetStandardElementsSol(void) {return standardElementsSol.data();}

inline unsigned short CMeshFEM_DG::GetNStandardMatchingFacesSol(void) {return standardMatchingFacesSol.size();}

inline FEMStandardInternalFaceClass* CMeshFEM_DG::GetStandardMatchingFacesSol(void) {return standardMatchingFacesSol.data();}
