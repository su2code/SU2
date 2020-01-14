/*!
 * \file fem_geometry_structure.inl
 * \brief In-Line subroutines of the <i>fem_geometry_structure.hpp</i> file.
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

inline CLong3T::CLong3T(){long0 = long1 = long2 = 0;}

inline CLong3T::CLong3T(const long a, const long b, const long c){long0 = a; long1 = b; long2 = c;}

inline CLong3T::~CLong3T(){}

inline CLong3T::CLong3T(const CLong3T &other){Copy(other);}

inline CLong3T& CLong3T::operator=(const CLong3T &other){Copy(other); return (*this);}

inline CReorderElements::~CReorderElements(void) { }

inline CReorderElements::CReorderElements(const CReorderElements &other) { Copy(other); }

inline CReorderElements& CReorderElements::operator=(const CReorderElements &other) { Copy(other); return (*this); }

inline bool CReorderElements::GetCommSolution(void) { return commSolution; }

inline unsigned short CReorderElements::GetElemType(void) { return elemType; }

inline unsigned long CReorderElements::GetGlobalElemID(void) { return globalElemID; }

inline unsigned short CReorderElements::GetTimeLevel(void) { return timeLevel; }

inline void CReorderElements::SetCommSolution(const bool val_CommSolution) { commSolution = val_CommSolution; }

inline CSortFaces::CSortFaces(unsigned long            val_nVolElemOwned,
                              unsigned long            val_nVolElemTot,
                              const CVolumeElementFEM *val_volElem)
 {nVolElemOwned = val_nVolElemOwned; nVolElemTot = val_nVolElemTot; volElem = val_volElem;}

inline CSortFaces::~CSortFaces(void) { }

inline CSortBoundaryFaces::CSortBoundaryFaces() { }

inline CSortBoundaryFaces::~CSortBoundaryFaces() { }

inline CVolumeElementFEM::CVolumeElementFEM(void) { }

inline CVolumeElementFEM::~CVolumeElementFEM(void) { }

inline CPointFEM::CPointFEM(void) { coor[0] = coor[1] = coor[2] = 0.0; }

inline CPointFEM::~CPointFEM(void) { }

inline CPointFEM::CPointFEM(const CPointFEM &other) { Copy(other); }

inline CPointFEM& CPointFEM::operator=(const CPointFEM &other) { Copy(other); return (*this); }

inline CInternalFaceElementFEM::CInternalFaceElementFEM(void) { }

inline CInternalFaceElementFEM::~CInternalFaceElementFEM(void) { }

inline CInternalFaceElementFEM::CInternalFaceElementFEM(const CInternalFaceElementFEM &other) { Copy(other); }

inline CInternalFaceElementFEM& CInternalFaceElementFEM::operator=(const CInternalFaceElementFEM &other) { Copy(other); return (*this); }

inline CSurfaceElementFEM::CSurfaceElementFEM(void) {indStandardElement = -1;}

inline CSurfaceElementFEM::~CSurfaceElementFEM(void) { }

inline CSurfaceElementFEM::CSurfaceElementFEM(const CSurfaceElementFEM &other) { Copy(other); }

inline CSurfaceElementFEM& CSurfaceElementFEM::operator=(const CSurfaceElementFEM &other) { Copy(other); return (*this); }

inline bool CSurfaceElementFEM::operator< (const CSurfaceElementFEM &other) const { return volElemID < other.volElemID; }

inline CBoundaryFEM::CBoundaryFEM(void) { periodicBoundary = haloInfoNeededForBC = false;  wallModel = NULL;}

inline CBoundaryFEM::~CBoundaryFEM(void) { if( wallModel ) delete wallModel; }

inline CMeshFEM::CMeshFEM(void) : CGeometry() { blasFunctions = NULL; }

inline CMeshFEM::~CMeshFEM(void) { if( blasFunctions ) {delete blasFunctions; blasFunctions = NULL;} }

inline CBoundaryFEM* CMeshFEM::GetBoundaries(void) {return boundaries.data();}

inline CPointFEM* CMeshFEM::GetMeshPoints(void) {return meshPoints.data();}

inline unsigned long CMeshFEM::GetNMeshPoints(void) {return meshPoints.size();}

inline unsigned long CMeshFEM::GetNVolElemOwned(void) {return nVolElemOwned;}

inline unsigned long CMeshFEM::GetNVolElemTot(void) {return nVolElemTot;}

inline CVolumeElementFEM* CMeshFEM::GetVolElem(void) {return volElem.data();}

inline unsigned long* CMeshFEM::GetNVolElemOwnedPerTimeLevel(void) {return nVolElemOwnedPerTimeLevel.data();}

inline unsigned long* CMeshFEM::GetNVolElemInternalPerTimeLevel(void) {return nVolElemInternalPerTimeLevel.data();}

inline unsigned long* CMeshFEM::GetNVolElemHaloPerTimeLevel(void) {return nVolElemHaloPerTimeLevel.data();}

inline vector<vector<unsigned long> > CMeshFEM::GetOwnedElemAdjLowTimeLevel(void) {return ownedElemAdjLowTimeLevel;}

inline vector<vector<unsigned long> > CMeshFEM::GetHaloElemAdjLowTimeLevel(void) {return haloElemAdjLowTimeLevel;}

inline unsigned short CMeshFEM::GetNStandardBoundaryFacesSol(void) {return standardBoundaryFacesSol.size();}

inline CFEMStandardBoundaryFace* CMeshFEM::GetStandardBoundaryFacesSol(void) {return standardBoundaryFacesSol.data();}

inline const vector<int>& CMeshFEM::GetRanksRecv(void) const {return ranksRecv;}

inline const vector<int>& CMeshFEM::GetRanksSend(void) const {return ranksSend;}

inline const vector<vector<unsigned long> >& CMeshFEM::GetEntitiesRecv(void) const {return entitiesRecv;}

inline const vector<vector<unsigned long> >& CMeshFEM::GetEntitiesSend(void) const {return entitiesSend;}

inline const vector<unsigned short>& CMeshFEM::GetRotPerMarkers(void) const {return rotPerMarkers;}

inline const vector<vector<unsigned long> >& CMeshFEM::GetRotPerHalos(void) const {return rotPerHalos;}

inline CMeshFEM_DG::CMeshFEM_DG(void) : CMeshFEM() { }

inline CMeshFEM_DG::~CMeshFEM_DG(void) { }

inline void CMeshFEM_DG::SetGlobal_to_Local_Point(void) {
  Global_to_Local_Point.clear();
  unsigned long ii = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {
    for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j, ++ii) {
      Global_to_Local_Point[volElem[i].offsetDOFsSolGlobal+j] = ii;
    }
  }
}

inline long CMeshFEM_DG::GetGlobal_to_Local_Point(unsigned long val_ipoint) const {
  auto it = Global_to_Local_Point.find(val_ipoint);
  if (it != Global_to_Local_Point.cend())
    return it->second;
  return -1;
}

inline su2double* CMeshFEM_DG::GetTimeCoefADER_DG(void) {return timeCoefADER_DG.data();}

inline su2double* CMeshFEM_DG::GetTimeInterpolDOFToIntegrationADER_DG(void) {return timeInterpolDOFToIntegrationADER_DG.data();}

inline su2double* CMeshFEM_DG::GetTimeInterpolAdjDOFToIntegrationADER_DG(void) {return timeInterpolAdjDOFToIntegrationADER_DG.data();}

inline unsigned long *CMeshFEM_DG::GetNMatchingFacesWithHaloElem(void) {return nMatchingFacesWithHaloElem.data();}

inline unsigned long *CMeshFEM_DG::GetNMatchingFacesInternal(void) {return nMatchingFacesInternal.data();}

inline CInternalFaceElementFEM* CMeshFEM_DG::GetMatchingFaces(void) {return matchingFaces.data();}

inline unsigned short CMeshFEM_DG::GetNStandardElementsSol(void) {return standardElementsSol.size();}

inline CFEMStandardElement* CMeshFEM_DG::GetStandardElementsSol(void) {return standardElementsSol.data();}

inline unsigned short CMeshFEM_DG::GetNStandardMatchingFacesSol(void) {return standardMatchingFacesSol.size();}

inline CFEMStandardInternalFace* CMeshFEM_DG::GetStandardMatchingFacesSol(void) {return standardMatchingFacesSol.data();}
