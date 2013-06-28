/*!
 * \file grid_movement_structure.inl
 * \brief In-Line subroutines of the <i>grid_movement_structure.hpp</i> file.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
 */
 
#pragma once

inline unsigned short CSurfaceMovement::GetnChunk(void) { return nChunk; }

inline bool CSurfaceMovement::GetChunkDefinition(void) { return ChunkDefinition; }

inline void CFreeFormChunk::SetnCornerPoints(unsigned short val_ncornerpoints){ nCornerPoints = val_ncornerpoints; }

inline unsigned short CFreeFormChunk::GetnCornerPoints(void){ return nCornerPoints; }

inline unsigned short CFreeFormChunk::GetnControlPoints(void){ return lOrder*mOrder*nOrder; }

inline unsigned long CFreeFormChunk::GetnSurfacePoints(void){ return 0; }

inline double *CFreeFormChunk::GetCoordCornerPoints(unsigned short val_icornerpoints) { return Coord_Corner_Points[val_icornerpoints]; }

inline double *CFreeFormChunk::GetCoordControlPoints(unsigned short val_iindex, unsigned short val_jindex, unsigned short val_kindex) { return Coord_Control_Points[val_iindex][val_jindex][val_kindex]; }

inline double  CFreeFormChunk::GetCoordCornerPoints(unsigned short val_dim, unsigned short val_icornerpoints) { return Coord_Corner_Points[val_icornerpoints][val_dim]; }

inline unsigned short CFreeFormChunk::GetlOrder(void) { return lOrder; }

inline unsigned short CFreeFormChunk::GetmOrder(void) { return mOrder; }

inline unsigned short CFreeFormChunk::GetnOrder(void) { return nOrder; }

inline void  CFreeFormChunk::SetCoordCornerPoints(double *val_coord, unsigned short val_icornerpoints) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++) 
		Coord_Corner_Points[val_icornerpoints][iDim] = val_coord[iDim];
}

inline void CFreeFormChunk::SetCoordControlPoints(double *val_coord, unsigned short iDegree, unsigned short jDegree, unsigned short kDegree) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++) 
		Coord_Control_Points[iDegree][jDegree][kDegree][iDim] = val_coord[iDim];
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

inline void CFreeFormChunk::CrossProduct (double *v1, double *v2, double *v3) {
	v3[0] = v1[1]*v2[2]-v1[2]*v2[1];
	v3[1] = v1[2]*v2[0]-v1[0]*v2[2];
	v3[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

inline double CFreeFormChunk::DotProduct (double *v1, double *v2) { double scalar = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]; return scalar; }

inline double CFreeFormChunk::GetNorm(double *a) { double  norm = sqrt(a[0]*a[0] + a[1]*a[1]+ a[2]*a[2]); return norm; }

inline void CFreeFormChunk::SetTag(string val_tag) { Tag = val_tag; }

inline string CFreeFormChunk::GetTag() { return Tag; }
