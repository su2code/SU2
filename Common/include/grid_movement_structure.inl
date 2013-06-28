/*!
 * \file grid_movement_structure.inl
 * \brief In-Line subroutines of the <i>grid_movement_structure.hpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.2
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

inline unsigned short CSurfaceMovement::GetnLevel(void) { return nLevel; }

inline unsigned short CSurfaceMovement::GetnChunk(void) { return nChunk; }

inline bool CSurfaceMovement::GetChunkDefinition(void) { return ChunkDefinition; }

inline void CFreeFormChunk::Set_MarkerIndex(unsigned short val_iMarker) { MarkerIndex.push_back(val_iMarker); }

inline void CFreeFormChunk::SetParentChunk(string val_iParentChunk) { ParentChunk.push_back(val_iParentChunk); }

inline void CFreeFormChunk::SetChildChunk(string val_iChildChunk) { ChildChunk.push_back(val_iChildChunk); }

inline void CFreeFormChunk::Set_VertexIndex(unsigned long val_iVertex) { VertexIndex.push_back(val_iVertex); }

inline void CFreeFormChunk::Set_PointIndex(unsigned long val_iPoint) { PointIndex.push_back(val_iPoint); }

inline void CFreeFormChunk::Set_CartesianCoord(double *val_coord) { CartesianCoord[0].push_back(val_coord[0]);
																																		CartesianCoord[1].push_back(val_coord[1]); 
																																		CartesianCoord[2].push_back(val_coord[2]); }
																																		
inline void CFreeFormChunk::Set_CartesianCoord(double *val_coord, unsigned long val_iSurfacePoints) { CartesianCoord[0][val_iSurfacePoints] = val_coord[0];
																																																			CartesianCoord[1][val_iSurfacePoints] = val_coord[1]; 
																																																			CartesianCoord[2][val_iSurfacePoints] = val_coord[2]; }		

inline void CFreeFormChunk::Set_ParametricCoord(double *val_coord) { ParametricCoord[0].push_back(val_coord[0]);
																																		 ParametricCoord[1].push_back(val_coord[1]); 
																																		 ParametricCoord[2].push_back(val_coord[2]); }
																																		 
inline void CFreeFormChunk::Set_ParametricCoord(double *val_coord, unsigned long val_iSurfacePoints) { ParametricCoord[0][val_iSurfacePoints] = val_coord[0];
																																																			 ParametricCoord[1][val_iSurfacePoints] = val_coord[1]; 
																																																			 ParametricCoord[2][val_iSurfacePoints] = val_coord[2]; }

inline unsigned short CFreeFormChunk::Get_MarkerIndex(unsigned long val_iSurfacePoints) { return MarkerIndex[val_iSurfacePoints]; }

inline unsigned short CFreeFormChunk::GetnParentChunk(void) { return ParentChunk.size(); }

inline unsigned short CFreeFormChunk::GetnChildChunk(void) { return ChildChunk.size(); }

inline string CFreeFormChunk::GetParentChunkTag(unsigned short val_ParentChunk) { return ParentChunk[val_ParentChunk]; }

inline string CFreeFormChunk::GetChildChunkTag(unsigned short val_ChildChunk) { return ChildChunk[val_ChildChunk]; }

inline unsigned long CFreeFormChunk::Get_VertexIndex(unsigned long val_iSurfacePoints) { return VertexIndex[val_iSurfacePoints]; }

inline unsigned long CFreeFormChunk::Get_PointIndex(unsigned long val_iSurfacePoints) { return PointIndex[val_iSurfacePoints]; }

inline double *CFreeFormChunk::Get_CartesianCoord(unsigned long val_iSurfacePoints) { 
																																										cart_coord_[0] = CartesianCoord[0][val_iSurfacePoints];
																																										cart_coord_[1] = CartesianCoord[1][val_iSurfacePoints];
																																										cart_coord_[2] = CartesianCoord[2][val_iSurfacePoints];
																																										return cart_coord_; }

inline double *CFreeFormChunk::Get_ParametricCoord(unsigned long val_iSurfacePoints) { 
																																										param_coord_[0] = ParametricCoord[0][val_iSurfacePoints];
																																										param_coord_[1] = ParametricCoord[1][val_iSurfacePoints];
																																										param_coord_[2] = ParametricCoord[2][val_iSurfacePoints];
																																										return param_coord_; }
																																										
inline unsigned long CFreeFormChunk::GetnSurfacePoint(void) { return PointIndex.size(); }

inline void CFreeFormChunk::SetnCornerPoints(unsigned short val_ncornerpoints){ nCornerPoints = val_ncornerpoints; }

inline unsigned short CFreeFormChunk::GetnCornerPoints(void){ return nCornerPoints; }

inline unsigned short CFreeFormChunk::GetnControlPoints(void){ return lOrder*mOrder*nOrder; }

inline unsigned long CFreeFormChunk::GetnSurfacePoints(void){ return 0; }

inline double *CFreeFormChunk::GetCoordCornerPoints(unsigned short val_icornerpoints) { return Coord_Corner_Points[val_icornerpoints]; }

inline double *CFreeFormChunk::GetCoordControlPoints(unsigned short val_iindex, unsigned short val_jindex, unsigned short val_kindex) { return Coord_Control_Points[val_iindex][val_jindex][val_kindex]; }

inline double *CFreeFormChunk::GetParCoordControlPoints(unsigned short val_iindex, unsigned short val_jindex, unsigned short val_kindex) { return ParCoord_Control_Points[val_iindex][val_jindex][val_kindex]; }

inline double  CFreeFormChunk::GetCoordCornerPoints(unsigned short val_dim, unsigned short val_icornerpoints) { return Coord_Corner_Points[val_icornerpoints][val_dim]; }

inline unsigned short CFreeFormChunk::GetlOrder(void) { return lOrder; }

inline unsigned short CFreeFormChunk::GetmOrder(void) { return mOrder; }

inline unsigned short CFreeFormChunk::GetnOrder(void) { return nOrder; }

inline void  CFreeFormChunk::SetCoordCornerPoints(double *val_coord, unsigned short val_icornerpoints) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++) 
		Coord_Corner_Points[val_icornerpoints][iDim] = val_coord[iDim];
}

inline void CFreeFormChunk::SetCoordControlPoints(double *val_coord, unsigned short iDegree, unsigned short jDegree, unsigned short kDegree) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++) {
			Coord_Control_Points[iDegree][jDegree][kDegree][iDim] = val_coord[iDim];
			Coord_Control_Points_Copy[iDegree][jDegree][kDegree][iDim] = Coord_Control_Points[iDegree][jDegree][kDegree][iDim];
		}
}

inline void CFreeFormChunk::SetParCoordControlPoints(double *val_coord, unsigned short iDegree, unsigned short jDegree, unsigned short kDegree) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
			ParCoord_Control_Points[iDegree][jDegree][kDegree][iDim] = val_coord[iDim];
}

inline void CFreeFormChunk::SetCoordCornerPoints(double val_xcoord, double val_ycoord, double val_zcoord, unsigned short val_icornerpoints) {
	Coord_Corner_Points[val_icornerpoints][0] = val_xcoord;
	Coord_Corner_Points[val_icornerpoints][1] = val_ycoord;
	Coord_Corner_Points[val_icornerpoints][2] = val_zcoord;
}

inline void CFreeFormChunk::SetControlPoints(unsigned short *val_index, double *movement) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Coord_Control_Points[val_index[0]][val_index[1]][val_index[2]][iDim] += movement[iDim];
}

inline void CFreeFormChunk::SetOriginalControlPoints() {
	for (unsigned short iDegree = 0; iDegree <= lDegree; iDegree++)
		for (unsigned short jDegree = 0; jDegree <= mDegree; jDegree++)
			for (unsigned short kDegree = 0; kDegree <= nDegree; kDegree++)
				for (unsigned short iDim = 0; iDim < nDim; iDim++)
					Coord_Control_Points[iDegree][jDegree][kDegree][iDim] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][iDim];
}

inline void CFreeFormChunk::CrossProduct (double *v1, double *v2, double *v3) {
	v3[0] = v1[1]*v2[2]-v1[2]*v2[1];
	v3[1] = v1[2]*v2[0]-v1[0]*v2[2];
	v3[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

inline double CFreeFormChunk::DotProduct (double *v1, double *v2) { double scalar = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]; return scalar; }

inline double CFreeFormChunk::GetNorm(double *a) { double  norm = sqrt(a[0]*a[0] + a[1]*a[1]+ a[2]*a[2]); return norm; }

inline void CFreeFormChunk::SetTag(string val_tag) { Tag = val_tag; }

inline string CFreeFormChunk::GetTag() { return Tag; }

inline void CFreeFormChunk::SetLevel(unsigned short val_level) { Level = val_level; }

inline unsigned short CFreeFormChunk::GetLevel() { return Level; }
