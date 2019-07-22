/*!
 * \file geometry_structure.inl
 * \brief In-Line subroutines of the <i>geometry_structure.hpp</i> file.
 * \author F. Palacios, T. Economon
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

inline CUnsignedLong2T::CUnsignedLong2T(){long0 = long1 = 0;}

inline CUnsignedLong2T::CUnsignedLong2T(const unsigned long a, const unsigned long b){long0 = a; long1 = b;}

inline CUnsignedLong2T::~CUnsignedLong2T(){}

inline CUnsignedLong2T::CUnsignedLong2T(const CUnsignedLong2T &other){Copy(other);}

inline CUnsignedLong2T& CUnsignedLong2T::operator=(const CUnsignedLong2T &other){Copy(other); return (*this);}

inline CUnsignedShort2T::CUnsignedShort2T(){short0 = short1 = 0;}

inline CUnsignedShort2T::CUnsignedShort2T(const unsigned short a, const unsigned short b){short0 = a; short1 = b;}

inline CUnsignedShort2T::~CUnsignedShort2T(){}

inline CUnsignedShort2T::CUnsignedShort2T(const CUnsignedShort2T &other){Copy(other);}

inline CUnsignedShort2T& CUnsignedShort2T::operator=(const CUnsignedShort2T &other){Copy(other); return (*this);}

inline CFaceOfElement::CFaceOfElement(const CFaceOfElement &other){Copy(other);}

inline CFaceOfElement& CFaceOfElement::operator=(const CFaceOfElement &other){Copy(other); return (*this);}

inline void CFaceOfElement::CreateUniqueNumbering(void){sort(cornerPoints, cornerPoints+nCornerPoints);}

inline CBoundaryFace::CBoundaryFace(const CBoundaryFace &other){Copy(other);}

inline CBoundaryFace& CBoundaryFace::operator=(const CBoundaryFace &other){Copy(other); return (*this);}

inline bool CBoundaryFace::operator<(const CBoundaryFace &other) const {return (globalBoundElemID < other.globalBoundElemID);}

inline CMatchingFace::CMatchingFace(const CMatchingFace &other){Copy(other);}

inline CMatchingFace& CMatchingFace::operator=(const CMatchingFace &other){Copy(other); return (*this);}

inline void CGeometry::SetGlobal_to_Local_Point(void) { }

inline long CGeometry::GetGlobal_to_Local_Point(unsigned long val_ipoint) { return 0; }

inline unsigned short CGeometry::GetGlobal_to_Local_Marker(unsigned short val_imarker) { return 0; }

inline unsigned long CGeometry::GetGlobal_nPoint(void) { return 0; }

inline void CGeometry::SetGlobal_nPointDomain(unsigned long val_global_npoint) { }

inline unsigned long CGeometry::GetGlobal_nPointDomain(void) { return 0; }

inline unsigned long CGeometry::GetGlobal_nElem(void) { return 0; }

inline unsigned long CGeometry::GetGlobal_nElemDomain(void) { return 0; }

inline unsigned long CGeometry::GetGlobal_nElemLine(void) { return 0; }

inline unsigned long CGeometry::GetGlobal_nElemTria(void) { return 0; }

inline unsigned long CGeometry::GetGlobal_nElemQuad(void) { return 0; }

inline unsigned long CGeometry::GetGlobal_nElemTetr(void) { return 0; }

inline unsigned long CGeometry::GetGlobal_nElemHexa(void) { return 0; }

inline unsigned long CGeometry::GetGlobal_nElemPris(void) { return 0; }

inline unsigned long CGeometry::GetGlobal_nElemPyra(void) { return 0; }

inline unsigned long CGeometry::GetnElemLine(void) { return 0; }

inline unsigned long CGeometry::GetnElemTria(void) { return 0; }

inline unsigned long CGeometry::GetnElemQuad(void) { return 0; }

inline unsigned long CGeometry::GetnElemTetr(void) { return 0; }

inline unsigned long CGeometry::GetnElemHexa(void) { return 0; }

inline unsigned long CGeometry::GetnElemPris(void) { return 0; }

inline unsigned long CGeometry::GetnElemPyra(void) { return 0; }

inline void CGeometry::Check_IntElem_Orientation(CConfig *config) { }

inline void CGeometry::Check_BoundElem_Orientation(CConfig *config) { }

inline void CGeometry::SetColorGrid(CConfig *config) { }

inline void CGeometry::SetColorGrid_Parallel(CConfig *config) { }

inline void CGeometry::SetColorFEMGrid_Parallel(CConfig *config) { }

inline void CGeometry::DivideConnectivity(CConfig *config, unsigned short Elem_Type) { }

inline void CGeometry::SetRotationalVelocity(CConfig *config, unsigned short val_iZone, bool print) { }

inline void CGeometry::SetShroudVelocity(CConfig *config){ }

inline void CGeometry::SetTranslationalVelocity(CConfig *config, unsigned short val_iZone, bool print) { }

inline void CGeometry::SetGridVelocity(CConfig *config, unsigned long iter) { }

inline void CGeometry::SetRestricted_GridVelocity(CGeometry *fine_mesh, CConfig *config) { }

inline void CGeometry::SetSendReceive(CConfig *config) { }

inline void CGeometry::SetBoundaries(CConfig *config) { }

inline void CGeometry::ComputeWall_Distance(CConfig *config) { }

inline void CGeometry::SetPositive_ZArea(CConfig *config) { }

inline void CGeometry::SetPoint_Connectivity(void) { }

inline void CGeometry::SetRCM_Ordering(CConfig *config) { }

inline void CGeometry::SetCoord_Smoothing (unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config) { }

inline void CGeometry::SetCoord(CGeometry *geometry) { }

inline void CGeometry::SetMultiGridWallHeatFlux(CGeometry *geometry, unsigned short val_marker){ }

inline void CGeometry::SetMultiGridWallTemperature(CGeometry *geometry, unsigned short val_marker){ }

inline void CGeometry::SetPoint_Connectivity(CGeometry *fine_grid) { }

inline void CGeometry::SetElement_Connectivity(void) { }

inline unsigned long CGeometry::GetnPoint(void) { return nPoint; }

inline unsigned long CGeometry::GetnPointDomain(void) { return nPointDomain; }

inline unsigned long CGeometry::GetnElem(void) { return nElem; }

inline unsigned short CGeometry::GetnDim(void) { return nDim; }

inline unsigned short CGeometry::GetnZone(void) { return nZone; }

inline unsigned short CGeometry::GetnMarker(void) { return nMarker; }

inline string CGeometry::GetMarker_Tag(unsigned short val_marker) { return Tag_to_Marker[val_marker]; }

inline unsigned long CGeometry::GetMax_GlobalPoint(void) { return Max_GlobalPoint; }

inline void CGeometry::SetnMarker(unsigned short val_nmarker) { nMarker = val_nmarker; }

inline void CGeometry::SetnElem_Bound(unsigned short val_marker, unsigned long val_nelem_bound) { nElem_Bound[val_marker]= val_nelem_bound; }

inline unsigned long CGeometry::GetnElem_Bound(unsigned short val_marker) { return nElem_Bound[val_marker]; }

inline void CGeometry::SetMarker_Tag(unsigned short val_marker, string val_index) { Tag_to_Marker[val_marker] = val_index; }

inline void CGeometry::SetnPoint(unsigned long val_npoint) { nPoint = val_npoint; }

inline void CGeometry::SetnPointDomain(unsigned long val_npoint) { nPointDomain = val_npoint; }

inline void CGeometry::SetnElem(unsigned long val_nelem) { nElem = val_nelem; }

inline void CGeometry::SetnDim(unsigned short val_nDim) { nDim = val_nDim; }

inline su2double* CGeometry::GetSpanWiseValue(unsigned short val_marker) { return SpanWiseValue[val_marker-1]; }

inline unsigned long CGeometry::GetnVertex(unsigned short val_marker) { return nVertex[val_marker]; }

inline unsigned short CGeometry::GetnSpanWiseSections(unsigned short marker_flag) { return nSpanWiseSections[marker_flag -1]; }

inline unsigned long CGeometry::GetnVertexSpan(unsigned short val_marker, unsigned short val_span) { return nVertexSpan[val_marker][val_span]; }

inline unsigned long CGeometry::GetnFreqSpan(unsigned short val_marker, unsigned short val_span) { return (nTotVertexSpan[val_marker][val_span]/2 -1); }

inline unsigned long CGeometry::GetnVertexSpanMax(unsigned short marker_flag){return nVertexSpanMax[marker_flag];}

inline unsigned long CGeometry::GetnFreqSpanMax(unsigned short marker_flag){return (nVertexSpanMax[marker_flag]/2 -1);}

inline void CGeometry::SetnVertexSpanMax(unsigned short marker_flag, unsigned long nVertMax){nVertexSpanMax[marker_flag] = nVertMax;}

inline unsigned long CGeometry::GetnEdge(void) { return nEdge; }

inline bool CGeometry::FindFace(unsigned long first_elem, unsigned long second_elem, unsigned short &face_first_elem, unsigned short &face_second_elem) { return 0;}

inline void CGeometry::SetBoundVolume(void) { }

inline void CGeometry::SetVertex(void) { }

inline void CGeometry::SetVertex(CConfig *config) { }

inline void CGeometry::ComputeNSpan(CConfig *config, unsigned short val_iZone, unsigned short marker_flag, bool allocate) { }

inline void CGeometry::SetTurboVertex(CConfig *config,unsigned short val_iZone, unsigned short marker_flag, bool allocate) { }

inline void CGeometry::UpdateTurboVertex(CConfig *config,unsigned short val_iZone, unsigned short marker_flag) { }

inline void CGeometry::SetAvgTurboValue(CConfig *config,unsigned short val_iZone, unsigned short marker_flag, bool allocate) { }

inline void CGeometry::GatherInOutAverageValues(CConfig *config, bool allocate){ }

inline su2double* CGeometry::GetAverageTurboNormal(unsigned short val_marker, unsigned short val_span){ return NULL;}

inline su2double* CGeometry::GetAverageGridVel(unsigned short val_marker, unsigned short val_span){ return NULL;}

inline su2double CGeometry::GetAverageTangGridVel(unsigned short val_marker, unsigned short val_span){return 0.0;}

inline su2double* CGeometry::GetAverageNormal(unsigned short val_marker, unsigned short val_span){return NULL;}

inline su2double CGeometry::GetSpanArea(unsigned short val_marker, unsigned short val_span){return 0.0;}

inline su2double CGeometry::GetTurboRadius(unsigned short val_marker, unsigned short val_span){return 0.0;}

inline su2double CGeometry::GetTangGridVelIn(unsigned short val_marker, unsigned short val_span){return 0.0;}

inline su2double CGeometry::GetTangGridVelOut(unsigned short val_marker, unsigned short val_span){return 0.0;}

inline su2double CGeometry::GetSpanAreaIn(unsigned short val_marker, unsigned short val_span){return 0.0;}

inline su2double CGeometry::GetSpanAreaOut(unsigned short val_marker, unsigned short val_span){return 0.0;}

inline su2double CGeometry::GetTurboRadiusIn(unsigned short val_marker, unsigned short val_span){return 0.0;}

inline su2double CGeometry::GetTurboRadiusOut(unsigned short val_marker, unsigned short val_span){return 0.0;}

inline void CGeometry::SetTangGridVelIn(su2double value, unsigned short val_marker, unsigned short val_span){ }

inline void CGeometry::SetTangGridVelOut(su2double value, unsigned short val_marker, unsigned short val_span){ }

inline void CGeometry::SetSpanAreaIn(su2double value, unsigned short val_marker, unsigned short val_span){ }

inline void CGeometry::SetSpanAreaOut(su2double value, unsigned short val_marker, unsigned short val_span){ }

inline void CGeometry::SetTurboRadiusIn(su2double value, unsigned short val_marker, unsigned short val_span){ }

inline void CGeometry::SetTurboRadiusOut(su2double value, unsigned short val_marker, unsigned short val_span){ }

inline unsigned long CGeometry::GetnTotVertexSpan(unsigned short val_marker, unsigned short val_span){return 0;}

inline su2double CGeometry::GetMinAngularCoord(unsigned short val_marker, unsigned short val_span){return 0.0;}

inline su2double CGeometry::GetMaxAngularCoord(unsigned short val_marker, unsigned short val_span){return 0.0;}

inline su2double CGeometry::GetMinRelAngularCoord(unsigned short val_marker, unsigned short val_span){return 0.0;}

inline void CGeometry::SetVertex(CGeometry *fine_grid, CConfig *config) { }

inline void CGeometry::SetCoord_CG(void) { }

inline void CGeometry::SetMaxLength(CConfig* config) { }

inline void CGeometry::SetControlVolume(CConfig *config, unsigned short action) { }

inline void CGeometry::SetControlVolume(CConfig *config, CGeometry *geometry, unsigned short action) { }

inline void CGeometry::VisualizeControlVolume(CConfig *config, unsigned short action) { }

inline void CGeometry::MatchNearField(CConfig *config) { }

inline void CGeometry::MatchActuator_Disk(CConfig *config) { }

inline void CGeometry::MatchInterface(CConfig *config) { }

inline void CGeometry::MatchPeriodic(CConfig *config, unsigned short val_periodic) { }

inline void CGeometry::SetBoundControlVolume(CConfig *config, unsigned short action) { }

inline void CGeometry::SetBoundControlVolume(CConfig *config, CGeometry *geometry, unsigned short action) { }

inline void CGeometry::SetTecPlot(char config_filename[MAX_STRING_SIZE], bool new_file) { }

inline void CGeometry::SetMeshFile(CConfig *config, string val_mesh_out_filename) { }

inline void CGeometry::SetMeshFile(CGeometry *geometry, CConfig *config, string val_mesh_out_filename) { }

inline void CGeometry::SetBoundTecPlot(char mesh_filename[MAX_STRING_SIZE], bool new_file, CConfig *config) { }

inline su2double CGeometry::Compute_MaxThickness(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) { return 0; }

inline su2double CGeometry::Compute_Twist(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) { return 0; }
  
inline su2double CGeometry::Compute_Chord(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) { return 0; }

inline su2double CGeometry::Compute_Width(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) { return 0; }

inline su2double CGeometry::Compute_WaterLineWidth(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) { return 0; }

inline su2double CGeometry::Compute_Height(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) { return 0; }

inline su2double CGeometry::Compute_LERadius(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) { return 0; }

inline su2double CGeometry::Compute_Thickness(su2double *Plane_P0, su2double *Plane_Normal, su2double Location, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil, su2double &ZLoc) { return 0; }

inline su2double CGeometry::Compute_Area(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) { return 0; }

inline su2double CGeometry::Compute_Length(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) { return 0; }

inline void CGeometry::Compute_Wing_LeadingTrailing(su2double *LeadingEdge, su2double *TrailingEdge, su2double *Plane_P0, su2double *Plane_Normal,
                                               vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil,
                                               vector<su2double> &Zcoord_Airfoil) { }

inline void CGeometry::Compute_Fuselage_LeadingTrailing(su2double *LeadingEdge, su2double *TrailingEdge, su2double *Plane_P0, su2double *Plane_Normal,
                                               vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil,
                                               vector<su2double> &Zcoord_Airfoil) { }

inline su2double CGeometry::Compute_Dihedral(su2double *LeadingEdge_im1, su2double *TrailingEdge_im1,
                                            su2double *LeadingEdge_i, su2double *TrailingEdge_i) { return 0; }

inline su2double CGeometry::Compute_Curvature(su2double *LeadingEdge_im1, su2double *TrailingEdge_im1,
                                              su2double *LeadingEdge_i, su2double *TrailingEdge_i,
                                              su2double *LeadingEdge_ip1, su2double *TrailingEdge_ip1) { return 0; }

inline void CGeometry::Compute_Wing(CConfig *config, bool original_surface,
                                    su2double &Wing_Volume, su2double &Wing_MinMaxThickness, su2double &Wing_MaxMaxThickness,
                                    su2double &Wing_MinChord, su2double &Wing_MaxChord,
                                    su2double &Wing_MinLERadius, su2double &Wing_MaxLERadius,
                                    su2double &Wing_MinToC, su2double &Wing_MaxToC, su2double &Wing_ObjFun_MinToC,
                                    su2double &Wing_MaxTwist, su2double &Wing_MaxCurvature,
                                    su2double &Wing_MaxDihedral) { }


inline void CGeometry::Compute_Fuselage(CConfig *config, bool original_surface,
  		                                su2double &Fuselage_Volume, su2double &Fuselage_WettedArea, 
  		                                su2double &Fuselage_MinWidth, su2double &Fuselage_MaxWidth,
  		                                su2double &Fuselage_MinWaterLineWidth, su2double &Fuselage_MaxWaterLineWidth,
  		                                su2double &Fuselage_MinHeight, su2double &Fuselage_MaxHeight,
  		                                su2double &Fuselage_MaxCurvature) { }

inline void CGeometry::Compute_Nacelle(CConfig *config, bool original_surface,
                                       su2double &Nacelle_Volume, su2double &Nacelle_MinMaxThickness, su2double &Nacelle_MaxMaxThickness,
                                       su2double &Nacelle_MinChord, su2double &Nacelle_MaxChord,
                                       su2double &Nacelle_MinLERadius, su2double &Nacelle_MaxLERadius,
                                       su2double &Nacelle_MinToC, su2double &Nacelle_MaxToC,
                                       su2double &Nacelle_ObjFun_MinToC, su2double &Nacelle_MaxTwist) { }

inline void CGeometry::FindNormal_Neighbor(CConfig *config) { }

inline void CGeometry::SetBoundSensitivity(CConfig *config) { }

inline su2double CGeometry::GetCustomBoundaryTemperature(unsigned short val_marker, unsigned long val_vertex){ return CustomBoundaryTemperature[val_marker][val_vertex]; }

inline void CGeometry::SetCustomBoundaryTemperature(unsigned short val_marker, unsigned long val_vertex, su2double val_customBoundaryTemperature){ CustomBoundaryTemperature[val_marker][val_vertex] = val_customBoundaryTemperature; }

inline su2double CGeometry::GetCustomBoundaryHeatFlux(unsigned short val_marker, unsigned long val_vertex){ return CustomBoundaryHeatFlux[val_marker][val_vertex]; }

inline void CGeometry::SetCustomBoundaryHeatFlux(unsigned short val_marker, unsigned long val_vertex, su2double val_customBoundaryHeatFlux){ CustomBoundaryHeatFlux[val_marker][val_vertex] = val_customBoundaryHeatFlux; }

inline void CGeometry::SetMGLevel(unsigned short val_iMesh) { MGLevel = val_iMesh; }

inline unsigned short CGeometry::GetMGLevel(void) { return MGLevel; }

inline void CPhysicalGeometry::SetPoint_Connectivity(CGeometry *geometry) { CGeometry::SetPoint_Connectivity(geometry); } 

inline void CMultiGridGeometry::SetPoint_Connectivity(void) { CGeometry::SetPoint_Connectivity(); }

inline void CPhysicalGeometry::SetGlobal_to_Local_Point(void) {
  Global_to_Local_Point.clear();
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Global_to_Local_Point[node[iPoint]->GetGlobalIndex()] = iPoint;
  }
}

inline long CPhysicalGeometry::GetGlobal_to_Local_Point(unsigned long val_ipoint) {
  map<unsigned long, unsigned long>::const_iterator MI = Global_to_Local_Point.find(val_ipoint);
  if (MI != Global_to_Local_Point.end()) {
    return Global_to_Local_Point[val_ipoint];
  } else {
    return -1;
  }
}

inline unsigned short CPhysicalGeometry::GetGlobal_to_Local_Marker(unsigned short val_imarker) { return Global_to_Local_Marker[val_imarker]; }

inline unsigned long CPhysicalGeometry::GetGlobal_nPoint(void) { return Global_nPoint; }

inline unsigned long CPhysicalGeometry::GetGlobal_nPointDomain(void) { return Global_nPointDomain; }

inline unsigned long CPhysicalGeometry::GetGlobal_nElem(void) { return Global_nElem; }

inline unsigned long CPhysicalGeometry::GetGlobal_nElemDomain(void) { return Global_nElemDomain; }

inline unsigned long CPhysicalGeometry::GetGlobal_nElemLine(void) { return Global_nelem_edge; }

inline unsigned long CPhysicalGeometry::GetGlobal_nElemTria(void) { return Global_nelem_triangle; }

inline unsigned long CPhysicalGeometry::GetGlobal_nElemQuad(void) { return Global_nelem_quad; }

inline unsigned long CPhysicalGeometry::GetGlobal_nElemTetr(void) { return Global_nelem_tetra; }

inline unsigned long CPhysicalGeometry::GetGlobal_nElemHexa(void) { return Global_nelem_hexa; }

inline unsigned long CPhysicalGeometry::GetGlobal_nElemPris(void) { return Global_nelem_prism; }

inline unsigned long CPhysicalGeometry::GetGlobal_nElemPyra(void) { return Global_nelem_pyramid; }

inline unsigned long CPhysicalGeometry::GetnElemLine(void) { return nelem_edge; }

inline unsigned long CPhysicalGeometry::GetnElemTria(void) { return nelem_triangle; }

inline unsigned long CPhysicalGeometry::GetnElemQuad(void) { return nelem_quad; }

inline unsigned long CPhysicalGeometry::GetnElemTetr(void) { return nelem_tetra; }

inline unsigned long CPhysicalGeometry::GetnElemHexa(void) { return nelem_hexa; }

inline unsigned long CPhysicalGeometry::GetnElemPris(void) { return nelem_prism; }

inline unsigned long CPhysicalGeometry::GetnElemPyra(void) { return nelem_pyramid; }

inline void CGeometry::SetGeometryPlanes(CConfig *config) {}

inline vector<su2double> CGeometry::GetGeometryPlanes() { return XCoordList; }

inline vector<su2double> CPhysicalGeometry::GetGeometryPlanes() { return XCoordList; }

inline vector<su2double> CMultiGridGeometry::GetGeometryPlanes() { return XCoordList; }

inline vector<vector<su2double> > CGeometry::GetXCoord() { return Xcoord_plane; }

inline vector<vector<su2double> > CPhysicalGeometry::GetXCoord() { return Xcoord_plane; }

inline vector<vector<su2double> > CMultiGridGeometry::GetXCoord() { return Xcoord_plane; }

inline vector<vector<su2double> > CGeometry::GetYCoord() { return Ycoord_plane; }

inline vector<vector<su2double> > CPhysicalGeometry::GetYCoord() { return Ycoord_plane; }

inline vector<vector<su2double> > CMultiGridGeometry::GetYCoord() { return Ycoord_plane; }

inline vector<vector<su2double> > CGeometry::GetZCoord() { return Zcoord_plane; }

inline vector<vector<su2double> > CPhysicalGeometry::GetZCoord() { return Zcoord_plane; }

inline vector<vector<su2double> > CMultiGridGeometry::GetZCoord() { return Zcoord_plane; }


inline vector<vector<unsigned long> > CGeometry::GetPlanarPoints() { return Plane_points; }

inline vector<vector<unsigned long> > CPhysicalGeometry::GetPlanarPoints() { return Plane_points; }

inline vector<vector<unsigned long> > CMultiGridGeometry::GetPlanarPoints() { return Plane_points; }

inline void CGeometry::SetSensitivity(CConfig* config) {}

inline void CGeometry::ReadUnorderedSensitivity(CConfig* config) {}

inline su2double CGeometry::GetSensitivity(unsigned long iPoint, unsigned short iDim) { return 0.0;}

inline su2double CPhysicalGeometry::GetSensitivity(unsigned long iPoint, unsigned short iDim) { return Sensitivity[iPoint*nDim+iDim];}

inline void CGeometry::SetSensitivity(unsigned long iPoint, unsigned short iDim, su2double val) {}

inline void CPhysicalGeometry::SetSensitivity(unsigned long iPoint, unsigned short iDim, su2double val) {Sensitivity[iPoint*nDim+iDim] = val;}

inline su2double* CPhysicalGeometry::GetAverageTurboNormal(unsigned short val_marker, unsigned short val_span){return AverageTurboNormal[val_marker][val_span];}

inline su2double* CPhysicalGeometry::GetAverageGridVel(unsigned short val_marker, unsigned short val_span){return AverageGridVel[val_marker][val_span]; }

inline su2double CPhysicalGeometry::GetAverageTangGridVel(unsigned short val_marker, unsigned short val_span){return AverageTangGridVel[val_marker][val_span]; }

inline su2double* CPhysicalGeometry::GetAverageNormal(unsigned short val_marker, unsigned short val_span){return AverageNormal[val_marker][val_span];}

inline su2double CPhysicalGeometry::GetSpanArea(unsigned short val_marker, unsigned short val_span){return SpanArea[val_marker][val_span];}

inline su2double CPhysicalGeometry::GetTurboRadius(unsigned short val_marker, unsigned short val_span){return TurboRadius[val_marker][val_span];}

inline su2double CPhysicalGeometry::GetTangGridVelIn(unsigned short val_marker, unsigned short val_span){return TangGridVelIn[val_marker][val_span];}

inline su2double CPhysicalGeometry::GetTangGridVelOut(unsigned short val_marker, unsigned short val_span){return TangGridVelOut[val_marker][val_span];}

inline su2double CPhysicalGeometry::GetSpanAreaIn(unsigned short val_marker, unsigned short val_span){return SpanAreaIn[val_marker][val_span];}

inline su2double CPhysicalGeometry::GetSpanAreaOut(unsigned short val_marker, unsigned short val_span){return SpanAreaOut[val_marker][val_span];}

inline su2double CPhysicalGeometry::GetTurboRadiusIn(unsigned short val_marker, unsigned short val_span){return TurboRadiusIn[val_marker][val_span];}

inline su2double CPhysicalGeometry::GetTurboRadiusOut(unsigned short val_marker, unsigned short val_span){return TurboRadiusOut[val_marker][val_span];}

inline void CPhysicalGeometry::SetTangGridVelIn(su2double value, unsigned short val_marker, unsigned short val_span){TangGridVelIn[val_marker][val_span] = value;}

inline void CPhysicalGeometry::SetTangGridVelOut(su2double value, unsigned short val_marker, unsigned short val_span){TangGridVelOut[val_marker][val_span] = value;}

inline void CPhysicalGeometry::SetSpanAreaIn(su2double value, unsigned short val_marker, unsigned short val_span){SpanAreaIn[val_marker][val_span] = value;}

inline void CPhysicalGeometry::SetSpanAreaOut(su2double value, unsigned short val_marker, unsigned short val_span){SpanAreaOut[val_marker][val_span] = value;}

inline void CPhysicalGeometry::SetTurboRadiusIn(su2double value, unsigned short val_marker, unsigned short val_span){ TurboRadiusIn[val_marker][val_span] = value;}

inline void CPhysicalGeometry::SetTurboRadiusOut(su2double value, unsigned short val_marker, unsigned short val_span){TurboRadiusOut[val_marker][val_span] = value;}

inline unsigned long CPhysicalGeometry::GetnTotVertexSpan(unsigned short val_marker, unsigned short val_span){return nTotVertexSpan[val_marker][val_span];}

inline su2double CPhysicalGeometry::GetMinAngularCoord(unsigned short val_marker, unsigned short val_span){return MinAngularCoord[val_marker][val_span];}

inline su2double CPhysicalGeometry::GetMaxAngularCoord(unsigned short val_marker, unsigned short val_span){return MaxAngularCoord[val_marker][val_span];}

inline su2double CPhysicalGeometry::GetMinRelAngularCoord(unsigned short val_marker, unsigned short val_span){return MinRelAngularCoord[val_marker][val_span];}

inline void CGeometry::Check_Periodicity(CConfig* config) {}
