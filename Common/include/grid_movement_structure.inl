/*!
 * \file grid_movement_structure.inl
 * \brief In-Line subroutines of the <i>grid_movement_structure.hpp</i> file.
 * \author F. Palacios, T. Economon, S. Padron
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

inline void CGridMovement::SetSurface_Deformation(CGeometry *geometry, CConfig *config) { }

inline unsigned short CSurfaceMovement::GetnLevel(void) { return nLevel; }

inline unsigned short CSurfaceMovement::GetnFFDBox(void) { return nFFDBox; }

inline bool CSurfaceMovement::GetFFDBoxDefinition(void) { return FFDBoxDefinition; }

inline void CFreeFormDefBox::Set_Fix_IPlane(unsigned short val_plane) { Fix_IPlane.push_back(val_plane); }

inline void CFreeFormDefBox::Set_Fix_JPlane(unsigned short val_plane) { Fix_JPlane.push_back(val_plane); }

inline void CFreeFormDefBox::Set_Fix_KPlane(unsigned short val_plane) { Fix_KPlane.push_back(val_plane); }

inline unsigned short CFreeFormDefBox::Get_Fix_IPlane(unsigned short val_index) { return Fix_IPlane[val_index]; }

inline unsigned short CFreeFormDefBox::Get_Fix_JPlane(unsigned short val_index) { return Fix_JPlane[val_index]; }

inline unsigned short CFreeFormDefBox::Get_Fix_KPlane(unsigned short val_index) { return Fix_KPlane[val_index]; }

inline unsigned short CFreeFormDefBox::Get_nFix_IPlane(void) { return Fix_IPlane.size(); }

inline unsigned short CFreeFormDefBox::Get_nFix_JPlane(void) { return Fix_JPlane.size(); }

inline unsigned short CFreeFormDefBox::Get_nFix_KPlane(void) { return Fix_KPlane.size(); }

inline void CFreeFormDefBox::Set_MarkerIndex(unsigned short val_iMarker) { MarkerIndex.push_back(val_iMarker); }

inline void CFreeFormDefBox::SetParentFFDBox(string val_iParentFFDBox) { ParentFFDBox.push_back(val_iParentFFDBox); }

inline void CFreeFormDefBox::SetChildFFDBox(string val_iChildFFDBox) { ChildFFDBox.push_back(val_iChildFFDBox); }

inline void CFreeFormDefBox::Set_VertexIndex(unsigned long val_iVertex) { VertexIndex.push_back(val_iVertex); }

inline void CFreeFormDefBox::Set_PointIndex(unsigned long val_iPoint) { PointIndex.push_back(val_iPoint); }

inline void CFreeFormDefBox::Set_CartesianCoord(su2double *val_coord) { CartesianCoord[0].push_back(val_coord[0]);
																																		CartesianCoord[1].push_back(val_coord[1]); 
																																		CartesianCoord[2].push_back(val_coord[2]); }
																																		
inline void CFreeFormDefBox::Set_CartesianCoord(su2double *val_coord, unsigned long val_iSurfacePoints) { CartesianCoord[0][val_iSurfacePoints] = val_coord[0];
																																																			CartesianCoord[1][val_iSurfacePoints] = val_coord[1]; 
																																																			CartesianCoord[2][val_iSurfacePoints] = val_coord[2]; }		

inline void CFreeFormDefBox::Set_ParametricCoord(su2double *val_coord) { ParametricCoord[0].push_back(val_coord[0]);
																																		 ParametricCoord[1].push_back(val_coord[1]); 
																																		 ParametricCoord[2].push_back(val_coord[2]); }
																																		 
inline void CFreeFormDefBox::Set_ParametricCoord(su2double *val_coord, unsigned long val_iSurfacePoints) { ParametricCoord[0][val_iSurfacePoints] = val_coord[0];
																																																			 ParametricCoord[1][val_iSurfacePoints] = val_coord[1]; 
																																																			 ParametricCoord[2][val_iSurfacePoints] = val_coord[2]; }

inline unsigned short CFreeFormDefBox::Get_MarkerIndex(unsigned long val_iSurfacePoints) { return MarkerIndex[val_iSurfacePoints]; }

inline unsigned short CFreeFormDefBox::GetnParentFFDBox(void) { return ParentFFDBox.size(); }

inline unsigned short CFreeFormDefBox::GetnChildFFDBox(void) { return ChildFFDBox.size(); }

inline string CFreeFormDefBox::GetParentFFDBoxTag(unsigned short val_ParentFFDBox) { return ParentFFDBox[val_ParentFFDBox]; }

inline string CFreeFormDefBox::GetChildFFDBoxTag(unsigned short val_ChildFFDBox) { return ChildFFDBox[val_ChildFFDBox]; }

inline unsigned long CFreeFormDefBox::Get_VertexIndex(unsigned long val_iSurfacePoints) { return VertexIndex[val_iSurfacePoints]; }

inline unsigned long CFreeFormDefBox::Get_PointIndex(unsigned long val_iSurfacePoints) { return PointIndex[val_iSurfacePoints]; }

inline su2double *CFreeFormDefBox::Get_CartesianCoord(unsigned long val_iSurfacePoints) { 
																																										cart_coord_[0] = CartesianCoord[0][val_iSurfacePoints];
																																										cart_coord_[1] = CartesianCoord[1][val_iSurfacePoints];
																																										cart_coord_[2] = CartesianCoord[2][val_iSurfacePoints];
																																										return cart_coord_; }

inline su2double *CFreeFormDefBox::Get_ParametricCoord(unsigned long val_iSurfacePoints) { 
																																										ParamCoord_[0] = ParametricCoord[0][val_iSurfacePoints];
																																										ParamCoord_[1] = ParametricCoord[1][val_iSurfacePoints];
																																										ParamCoord_[2] = ParametricCoord[2][val_iSurfacePoints];
																																										return ParamCoord_; }
																																										
inline unsigned long CFreeFormDefBox::GetnSurfacePoint(void) { return PointIndex.size(); }

inline void CFreeFormDefBox::SetnCornerPoints(unsigned short val_ncornerpoints) { nCornerPoints = val_ncornerpoints; }

inline unsigned short CFreeFormDefBox::GetnCornerPoints(void) { return nCornerPoints; }

inline unsigned short CFreeFormDefBox::GetnControlPoints(void) { return nControlPoints; }

inline void CFreeFormDefBox::SetnControlPoints(void) { nControlPoints = lOrder*mOrder*nOrder; }

inline unsigned long CFreeFormDefBox::GetnSurfacePoints(void) { return 0; }

inline su2double *CFreeFormDefBox::GetCoordCornerPoints(unsigned short val_icornerpoints) { return Coord_Corner_Points[val_icornerpoints]; }

inline su2double *CFreeFormDefBox::GetCoordControlPoints(unsigned short val_iindex, unsigned short val_jindex, unsigned short val_kindex) { return Coord_Control_Points[val_iindex][val_jindex][val_kindex]; }

inline su2double *CFreeFormDefBox::GetParCoordControlPoints(unsigned short val_iindex, unsigned short val_jindex, unsigned short val_kindex) { return ParCoord_Control_Points[val_iindex][val_jindex][val_kindex]; }

inline su2double  CFreeFormDefBox::GetCoordCornerPoints(unsigned short val_dim, unsigned short val_icornerpoints) { return Coord_Corner_Points[val_icornerpoints][val_dim]; }

inline unsigned short CFreeFormDefBox::GetlOrder(void) { return lOrder; }

inline unsigned short CFreeFormDefBox::GetmOrder(void) { return mOrder; }

inline unsigned short CFreeFormDefBox::GetnOrder(void) { return nOrder; }

inline void CFreeFormDefBox::SetlOrder(unsigned short val_lOrder) { lOrder = val_lOrder; lDegree = lOrder-1; }

inline void CFreeFormDefBox::SetmOrder(unsigned short val_mOrder) { mOrder = val_mOrder; mDegree = mOrder-1; }

inline void CFreeFormDefBox::SetnOrder(unsigned short val_nOrder) { nOrder = val_nOrder; nDegree = nOrder-1;}

inline void  CFreeFormDefBox::SetCoordCornerPoints(su2double *val_coord, unsigned short val_icornerpoints) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++) 
		Coord_Corner_Points[val_icornerpoints][iDim] = val_coord[iDim];
}

inline void CFreeFormDefBox::SetCoordControlPoints(su2double *val_coord, unsigned short iDegree, unsigned short jDegree, unsigned short kDegree) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++) {
			Coord_Control_Points[iDegree][jDegree][kDegree][iDim] = val_coord[iDim];
		}
}

inline void CFreeFormDefBox::SetCoordControlPoints_Copy(su2double *val_coord, unsigned short iDegree, unsigned short jDegree, unsigned short kDegree) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++) {
			Coord_Control_Points_Copy[iDegree][jDegree][kDegree][iDim] = val_coord[iDim];
		}
}

inline void CFreeFormDefBox::SetParCoordControlPoints(su2double *val_coord, unsigned short iDegree, unsigned short jDegree, unsigned short kDegree) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
			ParCoord_Control_Points[iDegree][jDegree][kDegree][iDim] = val_coord[iDim];
}

inline void CFreeFormDefBox::SetCoordCornerPoints(su2double val_xcoord, su2double val_ycoord, su2double val_zcoord, unsigned short val_icornerpoints) {
	Coord_Corner_Points[val_icornerpoints][0] = val_xcoord;
	Coord_Corner_Points[val_icornerpoints][1] = val_ycoord;
	Coord_Corner_Points[val_icornerpoints][2] = val_zcoord;
}

inline void CFreeFormDefBox::SetControlPoints(unsigned short *val_index, su2double *movement) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Coord_Control_Points[val_index[0]][val_index[1]][val_index[2]][iDim] += movement[iDim];
}

inline void CFreeFormDefBox::SetOriginalControlPoints() {
	for (unsigned short iDegree = 0; iDegree <= lDegree_Copy; iDegree++)
		for (unsigned short jDegree = 0; jDegree <= mDegree_Copy; jDegree++)
			for (unsigned short kDegree = 0; kDegree <= nDegree_Copy; kDegree++)
				for (unsigned short iDim = 0; iDim < nDim; iDim++)
					Coord_Control_Points[iDegree][jDegree][kDegree][iDim] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][iDim];
          
  lDegree = lDegree_Copy; mDegree = mDegree_Copy; nDegree = nDegree_Copy;
  lOrder = lOrder_Copy; mOrder = mOrder_Copy; nOrder = nOrder_Copy;
  nControlPoints = nControlPoints_Copy;
}

inline void CFreeFormDefBox::CrossProduct (su2double *v1, su2double *v2, su2double *v3) {
	v3[0] = v1[1]*v2[2]-v1[2]*v2[1];
	v3[1] = v1[2]*v2[0]-v1[0]*v2[2];
	v3[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

inline su2double CFreeFormDefBox::DotProduct (su2double *v1, su2double *v2) { su2double scalar = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]; return scalar; }

inline su2double CFreeFormDefBox::GetNorm(su2double *a) { su2double  norm = sqrt(a[0]*a[0] + a[1]*a[1]+ a[2]*a[2]); return norm; }

inline void CFreeFormDefBox::SetTag(string val_tag) { Tag = val_tag; }

inline string CFreeFormDefBox::GetTag() { return Tag; }

inline void CFreeFormDefBox::SetLevel(unsigned short val_level) { Level = val_level; }

inline unsigned short CFreeFormDefBox::GetLevel() { return Level; }

inline su2double CFreeFormDefBox::Determinant_3x3(su2double A00, su2double A01, su2double A02, su2double A10, su2double A11, su2double A12, su2double A20, su2double A21, su2double A22) {
	return A00*(A11*A22-A12*A21) - A01*(A10*A22-A12*A20) + A02*(A10*A21-A11*A20);
}

inline su2double CVolumetricMovement::Determinant_3x3(su2double A00, su2double A01, su2double A02, su2double A10, su2double A11, su2double A12, su2double A20, su2double A21, su2double A22) {
	return A00*(A11*A22-A12*A21) - A01*(A10*A22-A12*A20) + A02*(A10*A21-A11*A20);
}

inline void CVolumetricMovement::Set_nIterMesh(unsigned long val_nIterMesh) { nIterMesh = val_nIterMesh; }

inline unsigned long CVolumetricMovement::Get_nIterMesh() { return nIterMesh; }

inline void CVolumetricMovement::SetVolume_Deformation_Elas(CGeometry *geometry, CConfig *config, bool UpdateGeo, bool screen_output, bool Derivative) {  }

inline void CVolumetricMovement::Boundary_Dependencies(CGeometry **geometry, CConfig *config) {  }

inline bool CSurfaceMovement::CheckFFDBoxDefinition(CConfig *config, unsigned short iDV) {
  for (unsigned short iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
    if (FFDBox[iFFDBox]->GetTag() == config->GetFFDTag(iDV)) { return true;}
  }
  return false;
}

inline su2double CFreeFormBlending::GetBasis(short val_i, su2double val_t){return 0.0;}

inline su2double CFreeFormBlending::GetDerivative(short val_i, su2double val_t, short val_order){return 0.0;}

inline void CFreeFormBlending::SetOrder(short Order, short n_controlpoints){}

inline su2double CFreeFormBlending::GetOrder(){return Order;}

inline su2double CFreeFormBlending::GetDegree(){return Degree;}
