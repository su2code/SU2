/*!
 * \file interpolation_structure.hpp
 * \brief Headers of the main subroutines used by SU2_FSI.
 *        The subroutines and functions are in the <i>interpolation_structure.cpp</i> file.
 * \author H. Kline
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

#include "../../Common/include/mpi_structure.hpp"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>

#include "config_structure.hpp"
#include "geometry_structure.hpp"

using namespace std;


/*!
 * \class CInterpolator
 * \brief Main class for defining the interpolator, it requires
 * a child class for each particular interpolation method
 * \author H. Kline
 */
class CInterpolator {
protected:
  int rank, 	/*!< \brief MPI Rank. */
  size;       	/*!< \brief MPI Size. */
  unsigned int nZone;
  unsigned int donorZone, targetZone;

  unsigned long MaxLocalVertex_Donor,/*!\brief Maximum vertices per processor*/
  nGlobalFace_Donor,/*!\brief */
  nGlobalFaceNodes_Donor,/*!\brief */
  MaxFace_Donor,/*!\brief Maximum faces per processor*/
  MaxFaceNodes_Donor;/*!\brief Maximum nodes associated with faces per processor*/
  unsigned long *Buffer_Receive_nVertex_Donor, /*!\brief Buffer to store the number of vertices per processor on the Donor domain */
  *Buffer_Receive_nFace_Donor, /*!\brief Buffer to store the number of faces per processor*/
  *Buffer_Receive_nFaceNodes_Donor,/*!\brief Buffer to store the number of nodes associated with faces per processor*/
  *Buffer_Send_nVertex_Donor,/*!\brief Buffer to send number of vertices on the local processor*/
  *Buffer_Send_nFace_Donor,/*!\brief Buffer to send number of faces on the local processor*/
  *Buffer_Send_nFaceNodes_Donor,/*!\brief Buffer to send the number of nodes assocated with faces per processor*/
  *Buffer_Send_FaceIndex,/*!\brief Buffer to send indices pointing to the node indices that define the faces*/
  *Buffer_Receive_FaceIndex,/*!\brief Buffer to receive indices pointing to the node indices that define the faces*/
  *Buffer_Send_FaceNodes,/*!\brief Buffer to send indices pointing to the location of node information in other buffers, defining faces*/
  *Buffer_Receive_FaceNodes,/*!\brief Buffer to receive indices pointing to the location of node information in other buffers, defining faces*/
  *Buffer_Send_FaceProc,/*!\brief Buffer to send processor which stores the node indicated in Buffer_Receive_FaceNodes*/
  *Buffer_Receive_FaceProc;/*!\brief Buffer to receive processor which stores the node indicated in Buffer_Receive_FaceNodes*/

  long   *Buffer_Send_GlobalPoint,/*!\brief Buffer to send global point indices*/
  *Buffer_Receive_GlobalPoint; /*!\brief Buffer to receive global point indices*/

  su2double *Buffer_Send_Coord,/*!\brief Buffer to send coordinate values*/
  *Buffer_Send_Normal,/*!\brief Buffer to send normal vector values */
  *Buffer_Receive_Coord,/*!\brief Buffer to receive coordinate values*/
  *Buffer_Receive_Normal;/*!\brief Buffer to receive normal vector values*/
  
  unsigned long *Receive_GlobalPoint, /*!\brief Buffer to receive Global point indexes*/
  *Buffer_Receive_nLinkedNodes,       /*!\brief Buffer to receive the number of edges connected to each node*/
  *Buffer_Receive_LinkedNodes,        /*!\brief Buffer to receive the list of notes connected to the nodes through an edge*/
  *Buffer_Receive_StartLinkedNodes,   /*!\brief Buffer to receive the index of the Receive_LinkedNodes buffer where corresponding list of linked nodes begins */
  *Buffer_Receive_Proc;               /*!\brief Buffer to receive the thread that owns the node*/
  
  unsigned long  nGlobalVertex_Target, /*!\brief Global number of vertex of the target boundary*/
  nLocalVertex_Target,                 /*!\brief Number of vertex of the target boundary owned by the thread*/
  nGlobalVertex_Donor,                 /*!\brief Global number of vertex of the donor boundary*/
  nLocalVertex_Donor,                  /*!\brief Number of vertex of the donor boundary owned by the thread*/
  nGlobalVertex,                       /*!\brief Dummy variable to temporarily store the global number of vertex of a boundary*/
  nLocalLinkedNodes;                   /*!\brief Dummy variable to temporarily store the number of vertex of a boundary*/

public:
  CGeometry**** Geometry;        /*! \brief Vector which stores n zones of geometry. */
  CGeometry* donor_geometry;    /*! \brief Vector which stores the donor geometry. */
  CGeometry* target_geometry;   /*! \brief Vector which stores the target geometry. */

  /*!
   * \brief Constructor of the class.
   */
  CInterpolator(void);

  /*!
 * \brief Constructor of the class.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iZone - index of the donor zone
 * \param[in] jZone - index of the target zone
 */
  CInterpolator(CGeometry ****geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CInterpolator(void);

  /*!
   * \brief Find the index of the interface marker shared by that zone
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker_interface - Interface tag.
   */
  int Find_InterfaceMarker(CConfig *config, unsigned short val_marker_interface);

  /*!
   * \brief Check whether the interface should be processed or not
   * \param[in] val_markDonor  - Marker tag from donor zone.
   * \param[in] val_markTarget - Marker tag from target zone.
   */  
  bool CheckInterfaceBoundary(int val_markDonor, int val_markTarget);
  
  /*!
   * \brief Recontstruct the boundary connectivity from parallel partitioning and broadcasts it to all threads
   * \param[in] val_zone   - index of the zone
   * \param[in] val_marker - index of the marker
   */
  void ReconstructBoundary(unsigned long val_zone, int val_marker);
  
  /*!
   * \brief compute distance between 2 points
   * \param[in] point_i
   * \param[in] point_i
   */
  su2double PointsDistance(su2double *point_i, su2double *point_j);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Set_TransferCoeff(CConfig **config);

  /*!
   * \brief Determine array sizes used to collect and send coordinate and global point
   * information.
   * \param[in] faces - boolean that determines whether or not to set face information as well
   * \param[in] markDonor - Index of the boundary on the donor domain.
   * \param[in] markTarget - Index of the boundary on the target domain.
   * \param[in] nVertexDonor - Number of vertices on the donor boundary.
   * \param[in] nDim - number of physical dimensions.
   */
  void Determine_ArraySize(bool faces, int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim);

  /*!
   * \brief Collect and communicate vertex info: coord, global point, and if faces=true the normal vector
   * \param[in] faces - boolean that determines whether or not to set face information as well
   * \param[in] markDonor - Index of the boundary on the donor domain.
   * \param[in] markTarget - Index of the boundary on the target domain.
   * \param[in] nVertexDonor - Number of vertices on the donor boundary.
   * \param[in] nDim - number of physical dimensions.
   */
  void Collect_VertexInfo(bool faces, int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim);

};

/*!
 * \brief Nearest Neighbor interpolation
 */
class CNearestNeighbor : public CInterpolator {
public:

  /*!
   * \brief Constructor of the class.
   */
  CNearestNeighbor(void);

  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - index of the donor zone
   * \param[in] jZone - index of the target zone
   */
  CNearestNeighbor(CGeometry ****geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone);

  /*!
   * \brief Destructor of the class.
   */
  ~CNearestNeighbor(void);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   * \param[in] config - Definition of the particular problem.
   */
  void Set_TransferCoeff(CConfig **config);

};

/*!
 * \brief Isoparametric interpolation
  */
class CIsoparametric : public CInterpolator {
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - index of the donor zone
   * \param[in] jZone - index of the target zone
   */
  CIsoparametric(CGeometry ****geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone);

  /*!
   * \brief Destructor of the class.
   */
  ~CIsoparametric(void);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   * \param[in] config - Definition of the particular problem.
   */
  void Set_TransferCoeff(CConfig **config);

  /*!
   * \brief Calculate the isoparametric representation of point iVertex in marker iZone_0 by nodes of element donor_elem in marker jMarker of zone iZone_1.
   * \param[in] iVertex - vertex index of the point being interpolated.
   * \param[in] nDim - the dimension of the coordinates.
   * \param[in] iZone_1 - zone index of the element to use for interpolation (the DONOR zone)
   * \param[in] donor_elem - element index of the element to use for interpolation (or global index of a point in 2D)
   * \param[in[ nDonorPoints - number of donor points in the element.
   * \param[in] xj - point projected onto the plane of the donor element.
   * \param[out] isoparams - isoparametric coefficients. Must be allocated to size nNodes ahead of time. (size> nDonors)
   *
   * If the problem is 2D, the 'face' projected onto is actually an edge; the local index
   * of the edge is then stored in iFace, and the global index of the node (from which the edge
   * is referenced)
   */
  void Isoparameters(unsigned short nDim, unsigned short nDonor, su2double *X, su2double *xj,su2double* isoparams);

};

/*!
 * \brief Mirror interpolation: copy point linking and coefficient values from the opposing mesh
 * Assumes that the oppoosing mesh has already run interpolation. (otherwise this will result in empty/trivial interpolation)
 */
class CMirror : public CInterpolator {
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] geometry_container
   * \param[in] config - config container
   * \param[in] iZone - First zone
   * \param[in] jZone - Second zone
   *
   * Data is set in geometry[targetZone]
   *
   */
  CMirror(CGeometry ****geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone);

  /*!
   * \brief Destructor of the class.
   */
  ~CMirror(void);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   * \param[in] config - Definition of the particular problem.
   */
  void Set_TransferCoeff(CConfig **config);
  
};

/*!
 * \brief Sliding mesh approach
  */
class CSlidingMesh : public CInterpolator {
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - index of the donor zone
   * \param[in] jZone - index of the target zone
   */
  CSlidingMesh(CGeometry ****geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone);

  /*!
   * \brief Destructor of the class.
   */
  ~CSlidingMesh(void);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   * \param[in] config - Definition of the particular problem.
   */
  void Set_TransferCoeff(CConfig **config);
  
  /*!
   * \brief For 3-Dimensional grids, build the dual surface element
   * \param[in] map         - array containing the index of the boundary points connected to the node
   * \param[in] startIndex  - for each vertex specifies the corresponding index in the global array containing the indexes of all its neighbouring vertexes 
   * \param[in] nNeighbour  - for each vertex specifies the number of its neighbouring vertexes (on the boundary)
   * \param[in] coord       - array containing the coordinates of all the boundary vertexes
   * \param[in] centralNode - label of the vertex around which the dual surface element is built
   * \param[in] element  - double array where element node coordinates will be stored
   */  
  int Build_3D_surface_element(unsigned long *map, unsigned long *startIndex, unsigned long* nNeighbor, su2double *coord, unsigned long centralNode, su2double** element);
   
  /*!
   * \brief For 2-Dimensional grids, compute intersection length of two segments projected along a given direction
   * \param[in] A1 - first  point of segment A
   * \param[in] A2 - second point of segment A
   * \param[in] B1 - first  point of segment B
   * \param[in] B2 - second point of segment B
   * \param[in] Direction - along which segments are projected
   */
  su2double ComputeLineIntersectionLength(su2double* A1, su2double* A2, su2double* B1, su2double* B2, su2double* Direction);
  
  /*!
   * \brief For 3-Dimensional grids, compute intersection area between two triangle projected on a given plane
   * \param[in] A1 - first  point of triangle A
   * \param[in] A2 - second point of triangle A
   * \param[in] A3 - third  point of triangle A
   * \param[in] B1 - first  point of triangle B
   * \param[in] B2 - second point of triangle B
   * \param[in] B3 - third  point of triangle B
   * \param[in] Direction - vector normal to projection plane
   */
  su2double Compute_Triangle_Intersection(su2double* A1, su2double* A2, su2double* A3, su2double* B1, su2double* B2, su2double* B3, su2double* Direction);
  
  /*!
   * \brief For 3-Dimensional grids, compute intersection area between two triangle projected on a given plane
   * P1 from triangle P MUST be inside triangle Q, points order doesn't matter
   * \param[in] P1 - first  point of triangle A
   * \param[in] P2 - second point of triangle A
   * \param[in] P3 - third  point of triangle A
   * \param[in] Q1 - first  point of triangle B
   * \param[in] Q2 - second point of triangle B
   * \param[in] Q3 - third  point of triangle B
   */
  su2double ComputeIntersectionArea( su2double* P1, su2double* P2, su2double* P3, su2double* Q1, su2double* Q2, su2double* Q3 );
  
  /*!
   * \brief For 2-Dimensional grids, check whether, and compute, two lines are intersecting
   * \param[in] A1 - first  defining first line
   * \param[in] A2 - second defining first line
   * \param[in] B1 - first  defining second line
   * \param[in] B2 - second defining second line
   * \param[in] IntersectionPoint - Container for intersection coordinates
   */
  void ComputeLineIntersectionPoint( su2double* A1, su2double* A2, su2double* B1, su2double* B2, su2double* IntersectionPoint );
  
  /*!
   * \brief For N-Dimensional grids, check whether a point is inside a triangle specified by 3 T points
   * \param[in] Point - query point
   * \param[in] T1 - first  point of triangle T
   * \param[in] T2 - second point of triangle T
   * \param[in] T3 - third  point of triangle T
   */
  bool CheckPointInsideTriangle(su2double* Point, su2double* T1, su2double* T2, su2double* T3);
};

/*!
 * \brief Radial basis function interpolation
 */
class CRadialBasisFunction : public CInterpolator {
public:

  /*!
   * \brief Constructor of the class.
   */
  CRadialBasisFunction(void);

  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - index of the donor zone
   * \param[in] jZone - index of the target zone
   */
  CRadialBasisFunction(CGeometry ****geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone);

  /*!
   * \brief Destructor of the class.
   */
  ~CRadialBasisFunction(void);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   * \param[in] config - Definition of the particular problem.
   */
  void Set_TransferCoeff(CConfig **config);

  /*!
   * \brief Compute the value of a radial basis function, this is static so it can be re-used.
   * \param[in] type - of radial basis function
   * \param[in] radius - the characteristic dimension
   * \param[in] dist - distance
   */
  static su2double Get_RadialBasisValue(const short unsigned int type, const su2double &radius, const su2double &dist);
  
private:
  /*!
   * \brief If the polynomial term is included in the interpolation, and the points lie on a plane, the matrix becomes rank deficient
   * and cannot be inverted. This method detects that condition and corrects it by removing a row from P (the polynomial part of the matrix).
   * \param[in] m - number of rows of P
   * \param[in] n - number of columns of P
   * \param[in] skip_row - marks the row of P which is all ones (by construction)
   * \param[in] max_diff_tol_in - tolerance to detect points are on a plane
   * \param[out] keep_row - marks the rows of P kept
   * \param[out] n_polynomial - size of the polynomial part on exit (i.e. new number of rows)
   * \param[in,out] P - polynomial part of the matrix, may be changed or not!
   */
  void Check_PolynomialTerms(int m, unsigned long n, const int *skip_row, su2double max_diff_tol_in, int *keep_row, int &n_polynomial, su2double *P);

};

/*!
 * \brief Helper class used by CRadialBasisFunction to calculate the interpolation weights.
 * This does not inherit from CSysMatrix because: it is a dense format rather than block sparse;
 * as the interpolation is done on a single core there are no methods for communication.
 * The code can be compiled with LAPACK to use optimized matrix inversion and multiplication routines.
 * CPPFLAGS="-DHAVE_LAPACK" LDFLAGS=-L/path/to/lapack_lib LIBS="-llapack -lrefblas -lgfortran"
 */
class CSymmetricMatrix{

  private:
    
    bool initialized, inversed;
    int sz, num_val;
    int *perm_vec;
    passivedouble *val_vec, *decompose_vec, *inv_val_vec;

    enum DecompositionType { none, cholesky, lu };
    
    DecompositionType decomposed;
    
    inline int CalcIdx(int i, int j);
    inline int CalcIdxFull(int i, int j);
    inline void CheckBounds(int i, int j);
    
    passivedouble ReadL(int i, int j);
    passivedouble ReadU(int i, int j);
    passivedouble ReadInv(int i,int j);
    
    // not optimized dense matrix factorization and inversion for portability
    void CholeskyDecompose(bool overwrite);
    void LUDecompose();
    void CalcInv(bool overwrite);
    // matrix inversion using LAPACK routines (LDLT factorization)
    void CalcInv_sptri();
    void CalcInv_potri() {}; // LLT not implemented yet

  public:
	
    /*--- Methods ---*/
    CSymmetricMatrix();
    ~CSymmetricMatrix();

    void Initialize(int N);
    void Initialize(int N, su2double *formed_val_vec);

    inline int GetSize();

    void Write(int i, int j, const su2double& val);
    passivedouble Read(int i, int j);

    void MatVecMult(passivedouble *v);
    void MatMatMult(bool left_mult, su2double *mat_vec, int N);
    void Invert(const bool is_spd);

};
