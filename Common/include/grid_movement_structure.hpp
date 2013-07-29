/*!
 * \file grid_movement_structure.hpp
 * \brief Headers of the main subroutines for doing the numerical grid 
 *        movement (including volumetric movement, surface movement and Free From 
 *        technique definition). The subroutines and functions are in 
 *        the <i>grid_movement_structure.cpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
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

#include "geometry_structure.hpp"
#include "config_structure.hpp"
#include "matrix_structure.hpp"
#include "vector_structure.hpp"
#include "linear_solvers_structure.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;

/*! 
 * \class CGridMovement
 * \brief Class for moving the surface and volumetric 
 *        numerical grid (2D and 3D problems).
 * \author F. Palacios.
 * \version 2.0.6
 */
class CGridMovement {
public:
	
	/*! 
	 * \brief Constructor of the class. 
	 */
	CGridMovement(void);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CGridMovement(void);
  
  
  /*!
	 * \brief A pure virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetSurface_Deformation(CGeometry *geometry, CConfig *config);
  
};

/*! 
 * \class CFreeFormDefBox
 * \brief Class for defining the free form FFDBox structure.
 * \author F. Palacios & A. Galdran.
 * \version 2.0.6
 */
class CFreeFormDefBox : public CGridMovement {
public:
	unsigned short nDim;						/*!< \brief Number of dimensions of the problem. */
	unsigned short nCornerPoints,		/*!< \brief Number of corner points of the FFDBox. */
	nControlPoints;									/*!< \brief Number of control points of the FFDBox. */
	double **Coord_Corner_Points,		/*!< \brief Coordinates of the corner points. */
	****Coord_Control_Points,				/*!< \brief Coordinates of the control points. */
	****ParCoord_Control_Points,		/*!< \brief Coordinates of the control points. */
	****Coord_Control_Points_Copy,	/*!< \brief Coordinates of the control points (copy). */
	****Coord_SupportCP;						/*!< \brief Coordinates of the support control points. */
	unsigned short lOrder,	/*!< \brief Order of the FFDBox in the i direction. */
	mOrder,									/*!< \brief Order of the FFDBox in the j direction. */
	nOrder;									/*!< \brief Order of the FFDBox in the k direction. */
	unsigned short lDegree, /*!< \brief Degree of the FFDBox in the i direction. */
	mDegree,								/*!< \brief Degree of the FFDBox in the j direction. */
	nDegree;								/*!< \brief Degree of the FFDBox in the k direction. */
	double *param_coord, *param_coord_,	/*!< \brief Parametric coordinates of a point. */
	*cart_coord, *cart_coord_;			/*!< \brief Cartesian coordinates of a point. */
	double *gradient;			/*!< \brief Gradient of the point inversion process. */
	double MaxCoord[3];		/*!< \brief Maximum coordinates of the FFDBox. */
	double MinCoord[3];		/*!< \brief Minimum coordinates of the FFDBox. */
	string Tag;						/*!< \brief Tag to identify the FFDBox. */
	unsigned short Level;								/*!< \brief Nested level of the FFD box. */
	vector<double> CartesianCoord[3];		/*!< \brief Vector with all the cartesian coordinates in the FFD FFDBox. */
	vector<double> ParametricCoord[3];	/*!< \brief Vector with all the parametrics coordinates in the FFD FFDBox. */
	vector<unsigned short> MarkerIndex;	/*!< \brief Vector with all markers in the FFD FFDBox. */
	vector<unsigned long> VertexIndex;	/*!< \brief Vector with all vertex index in the FFD FFDBox. */
	vector<unsigned long> PointIndex;		/*!< \brief Vector with all points index in the FFD FFDBox. */
	unsigned long nSurfacePoint;				/*!< \brief Number of surfaces in the FFD FFDBox. */
	vector<string> ParentFFDBox;					/*!< \brief Vector with all the parent FFD FFDBox. */
	vector<string> ChildFFDBox;					/*!< \brief Vector with all the child FFD FFDBox. */
	
public:
	
	/*! 
	 * \brief Constructor of the class.
	 */
	CFreeFormDefBox(void);
	
	/*! 
	 * \overload
	 * \param[in] val_lDegree - Degree of the FFDBox in the i direction.
	 * \param[in] val_mDegree - Degree of the FFDBox in the j direction.
	 * \param[in] val_nDegree - Degree of the FFDBox in the k direction.
	 */	
	CFreeFormDefBox(unsigned short val_lDegree, unsigned short val_mDegree, unsigned short val_nDegree);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CFreeFormDefBox(void);
	
	/*! 
	 * \brief Add to the vector of markers a new marker.
	 * \param[in] val_iMarker - New marker inside the FFD box.
	 */	
	void Set_MarkerIndex(unsigned short val_iMarker);
	
	/*! 
	 * \brief Add to the vector of vertices a new vertex.
	 * \param[in] val_iVertex - New vertex inside the FFD box.
	 */	
	void Set_VertexIndex(unsigned long val_iVertex);
	
	/*! 
	 * \brief Add to the vector of points a new point.
	 * \param[in] val_iPoint - New point inside the FFD box.
	 */	
	void Set_PointIndex(unsigned long val_iPoint);
	
	/*! 
	 * \brief Add to the vector of cartesian coordinates a new coordinate.
	 * \param[in] val_coord - New coordinate inside the FFD box.
	 */	
	void Set_CartesianCoord(double *val_coord);
	
	/*! 
	 * \brief Add to the vector of parametric coordinates a new coordinate.
	 * \param[in] val_coord - New coordinate inside the FFD box.
	 */	
	void Set_ParametricCoord(double *val_coord);
	
	/*! 
	 * \brief Add to the vector of parent FFDBoxs a new FFD FFDBox.
	 * \param[in] val_iParentFFDBox - New parent FFDBox in the vector.
	 */	
	void SetParentFFDBox(string val_iParentFFDBox);
	
	/*! 
	 * \brief Add to the vector of child FFDBoxs a new FFD FFDBox.
	 * \param[in] val_iChildFFDBox - New child FFDBox in the vector.
	 */	
	void SetChildFFDBox(string val_iChildFFDBox);
	
	/*! 
	 * \brief _______________.
	 * \param[in] val_coord - _______________.
	 * \param[in] val_iSurfacePoints - _______________.
	 */	
	void Set_CartesianCoord(double *val_coord, unsigned long val_iSurfacePoints);
	
	/*! 
	 * \brief _______________.
	 * \param[in] val_coord - _______________.
	 * \param[in] val_iSurfacePoints - _______________.
	 */	
	void Set_ParametricCoord(double *val_coord, unsigned long val_iSurfacePoints);

	/*! 
	 * \brief _______________.
	 * \param[in] Get_MarkerIndex - _______________.
	 * \return _______________.
	 */	
	unsigned short Get_MarkerIndex(unsigned long val_iSurfacePoints);
	
	/*! 
	 * \brief _______________.
	 * \param[in] Get_VertexIndex - _______________.
	 * \return _______________.
	 */	
	unsigned long Get_VertexIndex(unsigned long val_iSurfacePoints);
	
	/*! 
	 * \brief _______________.
	 * \param[in] Get_PointIndex - _______________.
	 * \return _______________.
	 */	
	unsigned long Get_PointIndex(unsigned long val_iSurfacePoints);
	
	/*! 
	 * \brief _______________.
	 * \param[in] Get_CartesianCoord - _______________.
	 * \return _______________.
	 */	
	double *Get_CartesianCoord(unsigned long val_iSurfacePoints);
	
	/*! 
	 * \brief _______________.
	 * \param[in] Get_ParametricCoord - _______________.
	 * \return _______________.
	 */	
	double *Get_ParametricCoord(unsigned long val_iSurfacePoints);
	
	/*! 
	 * \brief _______________.
	 * \param[in] GetnSurfacePoint - _______________.
	 * \return _______________.
	 */	
	unsigned long GetnSurfacePoint(void);
	
	/*! 
	 * \brief _______________.
	 * \param[in] GetnParentFFDBox - _______________.
	 * \return _______________.
	 */	
	unsigned short GetnParentFFDBox(void);
	
	/*! 
	 * \brief _______________.
	 * \param[in] GetnChildFFDBox - _______________.
	 * \return _______________.
	 */	
	unsigned short GetnChildFFDBox(void);
	
	/*! 
	 * \brief _______________.
	 * \param[in] val_ParentFFDBox - _______________.
	 * \return _______________.
	 */	
	string GetParentFFDBoxTag(unsigned short val_ParentFFDBox);
	
	/*! 
	 * \brief _______________.
	 * \param[in] val_ChildFFDBox - _______________.
	 * \return _______________.
	 */	
	string GetChildFFDBoxTag(unsigned short val_ChildFFDBox);
	
	/*! 
	 * \brief Change the the position of the corners of the unitary FFDBox, 
	 *        and find the position of the control points for the FFDBox
	 * \param[in] FFDBox - Original FFDBox where we want to compute the control points.
	 */	
	void SetSupportCPChange(CFreeFormDefBox *FFDBox);
	
	/*! 
	 * \brief Set the number of corner points.
	 * \param[in] val_ncornerpoints - Number of corner points.
	 */	
	void SetnCornerPoints(unsigned short val_ncornerpoints);
	
	/*! 
	 * \brief Get the number of corner points.
	 * \return Number of corner points.
	 */	
	unsigned short GetnCornerPoints(void);
	
	/*! 
	 * \brief Get the number of control points.
	 * \return Number of control points.
	 */	
	unsigned short GetnControlPoints(void);
	
	/*! 
	 * \brief Get the number of numerical points on the surface.
	 * \return Number of numerical points on the surface.
	 */	
	unsigned long GetnSurfacePoints(void);
	
	/*! 
	 * \brief Set the corner point for the unitary FFDBox.
	 */	
	void SetUnitCornerPoints(void);
	
	/*! 
	 * \brief Set the coordinates of the corner points.
	 * \param[in] val_coord - Coordinates of the corner point with index <i>val_icornerpoints</i>.
	 * \param[in] val_icornerpoints - Index of the corner point.
	 */	
	void SetCoordCornerPoints(double *val_coord, unsigned short val_icornerpoints);

	/*! 
	 * \overload
	 * \param[in] val_xcoord - X coordinate of the corner point with index <i>val_icornerpoints</i>.
	 * \param[in] val_ycoord - Y coordinate of the corner point with index <i>val_icornerpoints</i>.
	 * \param[in] val_zcoord - Z coordinate of the corner point with index <i>val_icornerpoints</i>.
	 * \param[in] val_icornerpoints - Index of the corner point.
	 */	
	void SetCoordCornerPoints(double val_xcoord, double val_ycoord, double val_zcoord, unsigned short val_icornerpoints);
	
	/*! 
	 * \brief Set the coordinates of the control points.
	 * \param[in] val_coord - Coordinates of the control point.
	 * \param[in] iDegree - Index of the FFDBox, i direction.
	 * \param[in] jDegree - Index of the FFDBox, j direction.
	 * \param[in] kDegree - Index of the FFDBox, k direction.
	 */	
	void SetCoordControlPoints(double *val_coord, unsigned short iDegree, unsigned short jDegree, unsigned short kDegree);	

	/*! 
	 * \brief Set the coordinates of the control points.
	 * \param[in] val_coord - Coordinates of the control point.
	 * \param[in] iDegree - Index of the FFDBox, i direction.
	 * \param[in] jDegree - Index of the FFDBox, j direction.
	 * \param[in] kDegree - Index of the FFDBox, k direction.
	 */	
	void SetParCoordControlPoints(double *val_coord, unsigned short iDegree, unsigned short jDegree, unsigned short kDegree);	
	
	/*! 
	 * \brief Get the coordinates of the corner points.
	 * \param[in] val_dim - Index of the coordinate (x,y,z).
	 * \param[in] val_icornerpoints - Index of the corner point.
	 * \return Coordinate <i>val_dim</i> of the corner point <i>val_icornerpoints</i>.
	 */		
	double GetCoordCornerPoints(unsigned short val_dim, unsigned short val_icornerpoints);
	
	/*! 
	 * \brief Get the coordinates of the corner points.
	 * \param[in] val_icornerpoints - Index of the corner point.
	 * \return Pointer to the coordinate vector of the corner point <i>val_icornerpoints</i>.
	 */		
	double *GetCoordCornerPoints(unsigned short val_icornerpoints);
	
	/*! 
	 * \brief Get the coordinates of the control point.
	 * \param[in] val_iindex - Value of the local i index of the control point.
	 * \param[in] val_jindex - Value of the local j index of the control point.
	 * \param[in] val_kindex - Value of the local k index of the control point.
	 * \return Pointer to the coordinate vector of the control point with local index (i,j,k).
	 */		
	double *GetCoordControlPoints(unsigned short val_iindex, unsigned short val_jindex, unsigned short val_kindex);
	
	/*! 
	 * \brief Get the parametric coordinates of the control point.
	 * \param[in] val_iindex - Value of the local i index of the control point.
	 * \param[in] val_jindex - Value of the local j index of the control point.
	 * \param[in] val_kindex - Value of the local k index of the control point.
	 * \return Pointer to the coordinate vector of the control point with local index (i,j,k).
	 */		
	double *GetParCoordControlPoints(unsigned short val_iindex, unsigned short val_jindex, unsigned short val_kindex);
	
	/*! 
	 * \brief Set the control points in a parallelepiped (hexahedron).
	 */		
	void SetControlPoints_Parallelepiped(void);
	
	/*! 
	 * \brief Set the control points of the final chuck in a unitary hexahedron free form.
	 * \param[in] FFDBox - Original FFDBox where we want to compute the control points.
	 */	
	void SetSupportCP(CFreeFormDefBox *FFDBox);

	/*! 
	 * \brief Set the new value of the coordinates of the control points.
	 * \param[in] val_index - Local index (i,j,k) of the control point.
	 * \param[in] movement - Movement of the control point.
	 */	
	void SetControlPoints(unsigned short *val_index, double *movement);
	
	/*! 
	 * \brief Set the original value of the control points.
	 */	
	void SetOriginalControlPoints(void);
	
	/*! 
	 * \brief Set the tecplot file of the FFD chuck structure.
	 * \param[in] FFDBox_filename - Name of the output file with the FFD FFDBox structure.
	 * \param[in] new_file - New file or add to the existing file.
	 */		
	void SetTecplot(char FFDBox_filename[200], bool new_file);
	
	/*! 
	 * \brief Set the cartesian coords of a point in R^3 and convert them to the parametric coords of
	 *        our parametrization of a paralellepiped.
	 * \param[in] cart_coord - Cartesian coordinates of a point.
	 * \return Pointer to the parametric coordinates of a point.
	 */		
	double *GetParametricCoord_Analytical(double *cart_coord);
	
	/*! 
	 * \brief Iterative strategy for computing the parametric coordinates.
	 * \param[in] xyz - Cartesians coordinates of the target point.
	 * \param[in] guess - Initial guess for doing the parametric coordinates search.
	 * \param[in] tol - Level of convergence of the iterative method.
	 * \param[in] it_max - Maximal number of iterations.
	 * \return Parametric coordinates of the point.
	 * \attention Does this function work in boxes with degree greater than 9?. See ONERA M6 example
	 */		
	double *GetParametricCoord_Iterative(double *xyz, double *guess, double tol, unsigned long it_max);
	
	/*! 
	 * \brief Compute the cross product.
	 * \param[in] v1 - First input vector.
	 * \param[in] v2 - Second input vector.
	 * \param[out] v3 - Output vector wuth the cross product.
	 */		
	void CrossProduct(double *v1, double *v2, double *v3);
	
	/*! 
	 * \brief Compute the doc product.
	 * \param[in] v1 - First input vector.
	 * \param[in] v2 - Sencond input vector.
	 * \return Dot product between <i>v1</i>, and <i>v2</i>.
	 */		
	double DotProduct(double *v1, double *v2);
	
	/*! 
	 * \brief Here we take the parametric coords of a point in the box and we convert them to the 
	 *        physical cartesian coords by plugging the param_coords on the Bezier parameterization of our box.
	 * \param[in] param_coord - Parametric coordinates of a point.
	 * \return Pointer to the cartesian coordinates of a point.
	 */		
	double *EvalCartesianCoord(double *param_coord);
	
	/*! 
	 * \brief Set the Bernstein polynomial, defined as B_i^n(t) = Binomial(n,i)*t^i*(1-t)^(n-i).
	 * \param[in] val_n - Degree of the Bernstein polynomial.
	 * \param[in] val_i - Order of the Bernstein polynomial.
	 * \param[in] val_t - Value of the parameter where the polynomial is evaluated.
	 * \return Value of the Bernstein polynomial.
	 */		
	double GetBernstein(short val_n, short val_i, double val_t);
	
	/*! 
	 * \brief Get the binomial coefficient n over i, defined as n!/(m!(n-m)!)
	 * \note If the denominator is 0, the value is 1.
	 * \param[in] n - Upper coefficient.
	 * \param[in] m - Lower coefficient.
	 * \return Value of the binomial coefficient n over m.
	 */		
	unsigned short Binomial(unsigned short n, unsigned short m);
	
	/*! 
	 * \brief Get the binomial (optimized) coefficient n over m, defined as n!/(m!(n-m)!)
	 * \note If the denominator is 0, the value is 1.
	 * \param[in] n - Upper coefficient.
	 * \param[in] m - Lower coefficient.
	 * \return Value of the binomial coefficient n over m.
	 */		
	unsigned long BinomialOpt(unsigned long n, unsigned long m);

	/*! 
	 * \brief The Factorial Number n! is defined as n!=n*(n-1)*...*2*1.
	 * \param[in] n - Index of the factorial.
	 * \return Value of the factorial.
	 */		
	unsigned short Factorial(unsigned short n);
	
	/*! 
	 * \brief Get the order in the l direction of the FFD FFDBox.
	 * \return Order in the l direction of the FFD FFDBox.
	 */		
	unsigned short GetlOrder(void);
	
	/*! 
	 * \brief Get the order in the m direction of the FFD FFDBox.
	 * \return Order in the m direction of the FFD FFDBox.
	 */		
	unsigned short GetmOrder(void);
	
	/*! 
	 * \brief Get the order in the n direction of the FFD FFDBox.
	 * \return Order in the n direction of the FFD FFDBox.
	 */		
	unsigned short GetnOrder(void);
	
	/*! 
	 * \brief Set, at each vertex, the index of the free form FFDBox that contains the vertex.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iFFDBox - Index of the FFDBox.
	 */		
	bool GetPointFFD(CGeometry *geometry, CConfig *config, unsigned long iPoint);
	
	/*! 
	 * \brief Set the zone of the computational domain that is going to be deformed.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iFFDBox - Index of the FFDBox.
	 */		
	void SetDeformationZone(CGeometry *geometry, CConfig *config, unsigned short iFFDBox);
	
	/*! 
	 * \brief The "order" derivative of the i-th Bernstein polynomial of degree n, evaluated at t, 
	 *        is calculated as  (B_i^n(t))^{order}(t) = n*(GetBernstein(n-1,i-1,t)-GetBernstein(n-1,i,t)), 
	 *        having in account that if i=0, GetBernstein(n-1,-1,t) = 0.
	 * \param[in] val_n - Degree of the Bernstein polynomial.
	 * \param[in] val_i - Order of the Bernstein polynomial.
	 * \param[in] val_t - Value of the parameter where the polynomial is evaluated.
	 * \param[in] val_order - Order of the derivative.
	 * \return Value of the Derivative of the Bernstein polynomial.
	 */		
	double GetBernsteinDerivative(short val_n, short val_i, double val_t, short val_order);
	
	/*! 
	 * \brief The routine computes the gradient of F(u,v,w)=||X(u,v,w)-(x,y,z)||^2  evaluated at (u,v,w).
	 * \param[in] val_coord - Parametric coordiates of the target point.
	 * \param[in] xyz - Cartesians coordinates of the point.
	 * \return Value of the analytical gradient.
	 */		
	double *GetGradient_Analytical(double *val_coord, double *xyz);
	
	/*! 
	 * \brief The routine computes the numerical gradient of F(u,v,w)=||X(u,v,w)-(x,y,z)||^2  evaluated at (u,v,w).
	 * \param[in] uvw - Current value of the parametrics coordinates.
	 * \param[in] xyz - Cartesians coordinates of the target point to compose the functional.
	 * \return Value of the numerical gradient.
	 */		
	double *GetGradient_Numerical(double *uvw, double *xyz);
	
	/*! 
	 * \brief An auxiliary routine to help us compute the gradient of F(u,v,w)=||X(u,v,w)-(x,y,z)||^2 = 
	 *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)^2+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z)^2
	 *        Input: val_t, val_diff (to identify the index of the Bernstein polynomail we differentiate), the i,j,k , l,m,n 
	 *        E.G.: val_diff=2 => we differentiate w.r.t. w  (val_diff=0,1,or 2) Output: d [B_i^l*B_j^m *B_k^n] / d val_diff  
	 *        (val_u,val_v,val_w).
	 * \param[in] uvw - __________.
	 * \param[in] val_diff - __________.
	 * \param[in] ijk - __________.
	 * \param[in] lmn - Degree of the FFD box.
	 * \return __________.
	 */		
	double GetDerivative1(double *uvw, unsigned short val_diff, unsigned short *ijk, unsigned short *lmn);
	
	/*! 
	 * \brief An auxiliary routine to help us compute the gradient of F(u,v,w)=||X(u,v,w)-(x,y,z)||^2 =
	 *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)^2+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z)^2
	 *        Input: (u,v,w), dim , xyz=(x,y,z), l,m,n E.G.: dim=2 => we use the third coordinate of the control points, 
	 *        and the z-coordinate of xyz  (0<=dim<=2) Output: 2* ( (Sum_{i,j,k}^l,m,n P_{ijk}[dim] B_i^l[u] B_j^m[v] B_k^n[w]) - 
	 *        xyz[dim]).
	 * \param[in] uvw - __________.
	 * \param[in] dim - __________.
	 * \param[in] xyz - __________.
	 * \param[in] lmn - Degree of the FFD box.
	 * \return __________.
	 */		
	double GetDerivative2(double *uvw, unsigned short dim, double *xyz, unsigned short *lmn);
	
	/*! 
	 * \brief An auxiliary routine to help us compute the gradient of F(u,v,w)=||X(u,v,w)-(x,y,z)||^2 =
	 *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z)
	 * \param[in] uvw - Parametric coordiates of the point.
	 * \param[in] dim - Value of the coordinate to be differentiate.
	 * \param[in] diff_this - Diferentiation with respect this coordinate.
	 * \param[in] lmn - Degree of the FFD box.
	 * \return Sum_{i,j,k}^{l,m,n} [one of them with -1, 
	 *        depending on diff_this=0,1 or 2] P_{ijk}[dim] * (B_i^l[u] B_j^m[v] B_k^n[w])--one of them diffrentiated; 
	 *        which? diff_thiss will tell us ; E.G.: dim=2, diff_this=1 => we use the third coordinate of the control 
	 *        points, and derivate de v-Bersntein polynomial (use m-1 when summing!!).
	 */		
	double GetDerivative3(double *uvw, unsigned short dim, unsigned short diff_this,
						  unsigned short *lmn);
	
	/*! 
	 * \brief An auxiliary routine to help us compute the Hessian of F(u,v,w)=||X(u,v,w)-(x,y,z)||^2 =
	 *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z) 
	 *        Input: val_t, val_diff, val_diff2 (to identify the index of the Bernstein polynomials we differentiate), the i,j,k , l,m,n 
	 *        E.G.: val_diff=1, val_diff2=2  =>  we differentiate w.r.t. v and w  (val_diff=0,1,or 2)
	 *        E.G.: val_diff=0, val_diff2=0 => we differentiate w.r.t. u two times
	 *        Output: [d [B_i^l*B_j^m *B_k^n]/d val_diff *d [B_i^l*B_j^m *B_k^n]/d val_diff2] (val_u,val_v,val_w) .
	 * \param[in] uvw - __________.
	 * \param[in] val_diff - __________.
	 * \param[in] val_diff2 - __________.
	 * \param[in] ijk - __________.
	 * \param[in] lmn - Degree of the FFD box.
	 * \return __________.
	 */
	double GetDerivative4(double *uvw, unsigned short val_diff, unsigned short val_diff2,
						   unsigned short *ijk, unsigned short *lmn);
	
	/*! 
	 * \brief An auxiliary routine to help us compute the Hessian of F(u,v,w)=||X(u,v,w)-(x,y,z)||^2 =
	 *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z) 
	 *        Input: (u,v,w), dim , diff_this, diff_this_also, xyz=(x,y,z), l,m,n
	 *        Output:
	 *        Sum_{i,j,k}^{l,m,n} [two of them with -1, depending on diff_this,diff_this_also=0,1 or 2] 
	 *        P_{ijk}[dim] * (B_i^l[u] B_j^m[v] B_k^n[w])--one of them diffrentiated; which? diff_thiss will tell us ;
	 *        E.G.: dim=2, diff_this=1 => we use the third coordinate of the control points, and derivate de v-Bersntein 
	 *        polynomial (use m-1 when summing!!).
	 * \param[in] uvw - __________.
	 * \param[in] dim - __________.
	 * \param[in] diff_this - __________.
	 * \param[in] diff_this_also - __________.
	 * \param[in] lmn - Degree of the FFD box.
	 * \return __________.
	 */		
	double GetDerivative5(double *uvw, unsigned short dim, unsigned short diff_this, unsigned short diff_this_also,
						  unsigned short *lmn);
	
	/*! 
	 * \brief The routine that computes the Hessian of F(u,v,w)=||X(u,v,w)-(x,y,z)||^2 evaluated at (u,v,w)
	 *        Input: (u,v,w), (x,y,z)
	 *        Output: Hessian F (u,v,w).
	 * \param[in] uvw - Current value of the parametrics coordinates.
	 * \param[in] xyz - Cartesians coordinates of the target point to compose the functional.
	 * \param[in] val_Hessian - Value of the hessian.
	 */		
	void GetHessian_Analytical(double *uvw, double *xyz, double **val_Hessian);
	
	/*! 
	 * \brief Euclidean norm of a vector.
	 * \param[in] a - _______.
	 * \return __________.
	 */		
	double GetNorm(double *a);
	
	/*! 
	 * \brief Gauss method for solving a linear system.
	 * \param[in] A - __________.
	 * \param[in] rhs - __________.
	 * \param[in] nVar - __________.
	 */		
	void Gauss_Elimination(double** A, double* rhs, unsigned short nVar);
	
	/*! 
	 * \brief Set the tag that identify a FFDBox.
	 * \param[in] val_tag - value of the tag.
	 */	
	void SetTag(string val_tag);
	
	/*! 
	 * \brief Get the tag that identify a FFDBox.
	 * \return Value of the tag that identigy the FFDBox.
	 */	
	string GetTag(void);
	
	/*! 
	 * \brief Set the nested level of the FFDBox.
	 * \param[in] val_level - value of the level.
	 */	
	void SetLevel(unsigned short val_level);
	
	/*! 
	 * \brief Get the nested level of the FFDBox.
	 * \return Value of the nested level of the the FFDBox.
	 */	
	unsigned short GetLevel(void);
};

/*! 
 * \class CVolumetricMovement
 * \brief Class for moving the volumetric numerical grid.
 * \author F. Palacios, A. Bueno, T. Economon, S. Padron.
 * \version 2.0.6
 */
class CVolumetricMovement : public CGridMovement {
protected:

	unsigned short nDim;		/*!< \brief Number of dimensions. */
	unsigned short nVar;		/*!< \brief Number of variables. */
	unsigned long nPoint;		/*!< \brief Number of points. */
	unsigned long nPointDomain;		/*!< \brief Number of points in the domain. */

  CSysMatrix StiffMatrix; /*!< \brief Matrix to store the point-to-point stiffness. */
  CSysVector LinSysSol;
  CSysVector LinSysRes;

public:

	/*! 
	 * \brief Constructor of the class.
	 */
	CVolumetricMovement(CGeometry *geometry);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CVolumetricMovement(void);
	
	/*! 
	 * \brief __________________________
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void UpdateMultiGrid(CGeometry **geometry, CConfig *config);
	
	/*! 
	 * \brief Compute the stiffness matrix for grid deformation using spring analogy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \return Value of the length of the smallest edge of the grid.
	 */
	double SetSpringMethodContributions_Edges(CGeometry *geometry);
  
  /*!
	 * \brief Compute the stiffness matrix for grid deformation using spring analogy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \return Value of the length of the smallest edge of the grid.
	 */
	double SetFEAMethodContributions_Elem(CGeometry *geometry);
  
  /*!
	 * \brief Build the stiffness matrix for a 2-D triangular element. The result will be placed in StiffMatrix_Elem.
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] StiffMatrix_Elem - Element stiffness matrix to be filled.
	 * \param[in] val_Point_0 - Index value for Node 0 of the current triangle.
   * \param[in] val_Point_1 - Index value for Node 1 of the current triangle.
   * \param[in] val_Point_2 - Index value for Node 2 of the current triangle.
	 */
  void SetFEA_StiffMatrix2D(CGeometry *geometry, double **StiffMatrix_Elem,
                              unsigned long val_Point_0, unsigned long val_Point_1, unsigned long val_Point_2);
  
  /*!
	 * \brief Build the stiffness matrix for a 3-D tetrehedral element. The result will be placed in StiffMatrix_Elem.
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] StiffMatrix_Elem - Element stiffness matrix to be filled.
	 * \param[in] val_Point_0 - Index value for Node 0 of the current tetrahedron.
   * \param[in] val_Point_1 - Index value for Node 1 of the current tetrahedron.
   * \param[in] val_Point_2 - Index value for Node 2 of the current tetrahedron.
   * \param[in] val_Point_3 - Index value for Node 3 of the current tetrahedron.
	 */
  void SetFEA_StiffMatrix3D(CGeometry *geometry, double **StiffMatrix_Elem,
                            unsigned long val_Point_0, unsigned long val_Point_1, unsigned long val_Point_2, unsigned long val_Point_3);
	
  /*!
	 * \brief Add the stiffness matrix for a 2-D triangular element to the global stiffness matrix for the entire mesh (node-based).
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] StiffMatrix_Elem - Element stiffness matrix to be copied into the global matrix structure.
	 * \param[in] val_Point_0 - Index value for Node 0 of the current triangle.
   * \param[in] val_Point_1 - Index value for Node 1 of the current triangle.
   * \param[in] val_Point_2 - Index value for Node 2 of the current triangle.
	 */
  void AddFEA_StiffMatrix2D(CGeometry *geometry, double **StiffMatrix_Elem,
                              unsigned long val_Point_0, unsigned long val_Point_1, unsigned long val_Point_2);
  
  /*!
	 * \brief Add the stiffness matrix for a 2-D triangular element to the global stiffness matrix for the entire mesh (node-based).
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] StiffMatrix_Elem - Element stiffness matrix to be filled.
	 * \param[in] val_Point_0 - Index value for Node 0 of the current tetrahedron.
   * \param[in] val_Point_1 - Index value for Node 1 of the current tetrahedron.
   * \param[in] val_Point_2 - Index value for Node 2 of the current tetrahedron.
   * \param[in] val_Point_3 - Index value for Node 3 of the current tetrahedron.
	 */
  void AddFEA_StiffMatrix3D(CGeometry *geometry, double **StiffMatrix_Elem,
                            unsigned long val_Point_0, unsigned long val_Point_1, unsigned long val_Point_2, unsigned long val_Point_3);
  
  /*!
	 * \brief Check each of the mesh volumes for quality after performing grid deformation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */
	void CheckFEA_Grid(CGeometry *geometry);
  
  /*!
	 * \brief Check for negative volumes for 2-D elements after grid deformation.
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iElem - Index value for the current primal grid element.
	 * \param[in] val_Point_0 - Index value for Node 0 of the current triangle.
   * \param[in] val_Point_1 - Index value for Node 1 of the current triangle.
   * \param[in] val_Point_2 - Index value for Node 2 of the current triangle.
	 */
  void CheckFEA_Elem2D(CGeometry *geometry, unsigned long val_iElem, unsigned long val_Point_0, unsigned long val_Point_1, unsigned long val_Point_2);
  
  /*!
	 * \brief Check for negative volumes for 3-D elements after grid deformation.
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iElem - Index value for the current primal grid element.
	 * \param[in] val_Point_0 - Index value for Node 0 of the current tetrahedron.
   * \param[in] val_Point_1 - Index value for Node 1 of the current tetrahedron.
   * \param[in] val_Point_2 - Index value for Node 2 of the current tetrahedron.
   * \param[in] val_Point_3 - Index value for Node 3 of the current tetrahedron.
	 */
  void CheckFEA_Elem3D(CGeometry *geometry, unsigned long val_iElem, unsigned long val_Point_0, unsigned long val_Point_1, unsigned long val_Point_2, unsigned long val_Point_3);
  
	/*!
	 * \brief Check the boundary vertex that are going to be moved.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetBoundaryDisplacements(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Check the domain points vertex that are going to be moved.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetDomainDisplacements(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Update the value of the coordinates after the the grid movement.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void UpdateGridCoord(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Unsteady grid movement using rigid mesh rotation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Physical time iteration number.
	 */
	void SetRigidRotation(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter);
	
  /*!
	 * \brief Unsteady pitching grid movement using rigid mesh motion.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Physical time iteration number.
	 */
	void SetRigidPitching(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter);
  
  /*!
	 * \brief Unsteady plunging grid movement using rigid mesh motion.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Physical time iteration number.
	 */
	void SetRigidPlunging(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter);
  
  /*!
	 * \brief Unsteady translational grid movement using rigid mesh motion.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Physical time iteration number.
	 */
	void SetRigidTranslation(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter);

  /*!
	 * \brief Solve the typical section wing model.
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Cl - Coefficient of lift at particular iteration.
   * \param[in] Cm - Moment coefficient about z-axis at particular iteration.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Physical time iteration number.
   * \param[in] displacements - solution of typical section wing model.
	 */
  void SolveTypicalSectionWingModel(CGeometry *geometry, double Cl, double Cm, CConfig *config, unsigned short iZone, unsigned long iter, double (&displacements)[4]);
  
    /*!
	 * \brief Unsteady aeroelastic grid movement by deforming the mesh.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
     * \param[in] iZone - Zone number in the mesh.
     * \param[in] displacements - solution of typical section wing model.
	 */
    void AeroelasticDeform(CGeometry *geometry, CConfig *config, unsigned short iZone, double displacements[4]);
    
    /*!
	 * \brief Unsteady aeroelastic grid movement using rigid mesh motion.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
     * \param[in] iZone - Zone number in the mesh.
     * \param[in] displacements - solution of typical section wing model.
	 */
  void AeroelasticRigid(CGeometry *geometry, CConfig *config, unsigned short iZone, double displacements[4]);
  
  /*!
	 * \brief Sets up the generalized eigenvectors and eigenvalues needed to solve the structural equations.
   * \param[in] PHI - Matrix of the generalized eigenvectors.
   * \param[in] lambda - The eigenvalues of the generalized eigensystem.
	 * \param[in] config - Definition of the particular problem.
	 */
  void SetUpTypicalSectionWingModel(double (&PHI)[2][2],double (&lambda)[2], CConfig *config);
  
  /*!
	 * \brief Grid deformation using the spring analogy method.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] UpdateGeo - Update geometry.
	 */
	void SetVolume_Deformation(CGeometry *geometry, CConfig *config, bool UpdateGeo);
  
  /*!
	 * \brief Compute the determinant of a 3 by 3 matrix.
	 * \param[in] val_matrix 3 by 3 matrix.
	 * \result Determinant of the matrix
	 */
	double Determinant_3x3(double A00, double A01, double A02, double A10, double A11, double A12, double A20, double A21, double A22);

};

/*! 
 * \class CSurfaceMovement
 * \brief Class for moving the surface numerical grid.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CSurfaceMovement : public CGridMovement {
protected:
  CFreeFormDefBox** FFDBox;	/*!< \brief Definition of the Free Form Deformation Box. */
	unsigned short nFFDBox;	/*!< \brief Number of FFD FFDBoxs. */
	unsigned short nLevel;	/*!< \brief Level of the FFD FFDBoxs (parent/child). */
	bool FFDBoxDefinition;	/*!< \brief If the FFD FFDBox has been defined in the input file. */

public:
	
	/*! 
	 * \brief Constructor of the class.
	 */
	CSurfaceMovement(void);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSurfaceMovement(void);
	
	/*! 
	 * \brief Set a Hicks-Henne deformation bump functions on an airfoil.
	 * \param[in] boundary - Geometry of the boundary.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iDV - Index of the design variable.
	 * \param[in] ResetDef - Reset the deformation before starting a new one.
	 */
	void SetHicksHenne(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef);
  
  /*!
	 * \brief Set a spherical design problem.
	 * \param[in] boundary - Geometry of the boundary.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iDV - Index of the design variable.
	 * \param[in] ResetDef - Reset the deformation before starting a new one.
	 */
	void SetSpherical(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef);
	
  /*!
	 * \brief Set a Hicks-Henne deformation bump functions on an airfoil.
	 * \param[in] boundary - Geometry of the boundary.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iDV - Index of the design variable.
	 * \param[in] ResetDef - Reset the deformation before starting a new one.
	 */
	void SetCosBump(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef);
  
  /*!
	 * \brief Set a Hicks-Henne deformation bump functions on an airfoil.
	 * \param[in] boundary - Geometry of the boundary.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iDV - Index of the design variable.
	 * \param[in] ResetDef - Reset the deformation before starting a new one.
	 */
	void SetFourier(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef);
  
	/*! 
	 * \brief Set a NACA 4 digits airfoil family for airfoil deformation.
	 * \param[in] boundary - Geometry of the boundary.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetNACA_4Digits(CGeometry *boundary, CConfig *config);
	
	/*! 
	 * \brief Set a parabolic family for airfoil deformation.
	 * \param[in] boundary - Geometry of the boundary.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetParabolic(CGeometry *boundary, CConfig *config);
	
	/*! 
	 * \brief Set a obstacle in a channel.
	 * \param[in] boundary - Geometry of the boundary.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetObstacle(CGeometry *boundary, CConfig *config);
	
	/*! 
	 * \brief Stretch one side of a channel.
	 * \param[in] boundary - Geometry of the boundary.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetStretch(CGeometry *boundary, CConfig *config);
	
	/*! 
	 * \brief Set a rotation for surface movement.
	 * \param[in] boundary - Geometry of the boundary.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iDV - Index of the design variable.
	 * \param[in] ResetDef - Reset the deformation before starting a new one.
	 */
	void SetRotation(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef);
  
  /*! 
	 * \brief Deforms a 2-D flutter/pitching surface during an unsteady simulation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iter - Current physical time iteration.
	 */
	void SetBoundary_Flutter2D(CGeometry *geometry, CConfig *config, 
                             unsigned long iter);
	
	/*! 
	 * \brief Deforms a 3-D flutter/pitching surface during an unsteady simulation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iter - Current physical time iteration.
	 */
	void SetBoundary_Flutter3D(CGeometry *geometry, CConfig *config, 
                             CFreeFormDefBox **FFDBox, unsigned long iter);	
	
  /*! 
	 * \brief Set the collective pitch for a blade surface movement.
   * \author T. Economon
	 * \param[in] boundary - Geometry of the boundary.
	 * \param[in] config - Definition of the particular problem.
	 */
  void SetCollective_Pitch(CGeometry *geometry, CConfig *config);
  
  /*! 
	 * \brief Set any surface deformationsbased on an input file.
   * \author T. Economon
	 * \param[in] boundary - Geometry of the boundary.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Current physical time iteration.
	 */
  void SetExternal_Deformation(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter);
  
	/*! 
	 * \brief Set a displacement for surface movement.
	 * \param[in] boundary - Geometry of the boundary.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iDV - Index of the design variable.
	 * \param[in] ResetDef - Reset the deformation before starting a new one.
	 */
	void SetDisplacement(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef);
	
	/*! 
	 * \brief Copy the boundary coordinates to each vertex.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void CopyBoundary(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Set the surface/boundary deformation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetSurface_Deformation(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Compute the parametric coordinates of a grid point using a point inversion strategy
	 *        in the free form FFDBox.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] FFDBox - Array with all the free forms FFDBoxs of the computation.
	 */		
	void SetParametricCoord(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox);
	
	/*! 
	 * \brief Update the parametric coordinates of a grid point using a point inversion strategy
	 *        in the free form FFDBox.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] FFDBox - Array with all the free forms FFDBoxs of the computation.
	 */		
	void UpdateParametricCoord(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox);
	
	/*! 
	 * \brief _____________________.
	 * \param[in] geometry - _____________________.
	 * \param[in] config - _____________________.
	 * \param[in] FFDBox - _____________________.
	 */	
	void SetParametricCoordCP(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBoxParent, CFreeFormDefBox *FFDBoxChild);
	
	/*! 
	 * \brief _____________________.
	 * \param[in] geometry - _____________________.
	 * \param[in] config - _____________________.
	 * \param[in] FFDBox - _____________________.
	 */	
	void GetCartesianCoordCP(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBoxParent, CFreeFormDefBox *FFDBoxChild);

	/*! 
	 * \brief Recompute the cartesian coordinates using the control points position.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] FFDBox - Array with all the free forms FFDBoxs of the computation.
	 */		
	void SetCartesianCoord(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox);
	
	/*! 
	 * \brief Set the deformation of the Free From box using the control point position.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] FFDBox - Array with all the free forms FFDBoxs of the computation.
	 * \param[in] iDV - Index of the design variable.
	 * \param[in] ResetDef - Reset the deformation before starting a new one.
	 */		
	void SetFFDCPChange(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, unsigned short iDV, bool ResetDef);
	
	/*! 
	 * \brief Set a camber deformation of the Free From box using the control point position.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] FFDBox - Array with all the free forms FFDBoxs of the computation.
	 * \param[in] iDV - Index of the design variable.
	 * \param[in] ResetDef - Reset the deformation before starting a new one.
	 */		
	void SetFFDCamber(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, unsigned short iDV, bool ResetDef);
	
	/*! 
	 * \brief Set a thickness deformation of the Free From box using the control point position.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] FFDBox - Array with all the free forms FFDBoxs of the computation.
	 * \param[in] iDV - Index of the design variable.
	 * \param[in] ResetDef - Reset the deformation before starting a new one.
	 */		
	void SetFFDThickness(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, unsigned short iDV, bool ResetDef);
	
	/*! 
	 * \brief Set a volume deformation of the Free From box using the control point position.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] FFDBox - Array with all the free forms FFDBoxs of the computation.
	 * \param[in] iDV - Index of the design variable.
	 * \param[in] ResetDef - Reset the deformation before starting a new one.
	 */		
	void SetFFDVolume(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, unsigned short iDV, bool ResetDef);
	
	/*! 
	 * \brief Set a dihedral angle deformation of the Free From box using the control point position.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] FFDBox - Array with all the free forms FFDBoxs of the computation.
	 * \param[in] iDV - Index of the design variable.
	 * \param[in] ResetDef - Reset the deformation before starting a new one.
	 */		
	void SetFFDDihedralAngle(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, unsigned short iDV, bool ResetDef);
	
	/*! 
	 * \brief Set a twist angle deformation of the Free From box using the control point position.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] FFDBox - Array with all the free forms FFDBoxs of the computation.
	 * \param[in] iDV - Index of the design variable.
	 * \param[in] ResetDef - Reset the deformation before starting a new one.
	 */		
	void SetFFDTwistAngle(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, unsigned short iDV, bool ResetDef);
	
	/*! 
	 * \brief Set a rotation angle deformation of the Free From box using the control point position.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] FFDBox - Array with all the free forms FFDBoxs of the computation.
	 * \param[in] iDV - Index of the design variable.
	 * \param[in] ResetDef - Reset the deformation before starting a new one.
	 */		
	void SetFFDRotation(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, unsigned short iDV, bool ResetDef);
	
	/*! 
	 * \brief Read the free form information from the grid input file.
	 * \note If there is no control point information, and no parametric 
	 *       coordinates information, the code will compute that information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] FFDBox - Array with all the free forms FFDBoxs of the computation.
	 * \param[in] val_mesh_filename - Name of the grid input file.
   * \param[in] val_vertex - With vertex information.
	 */		
	void ReadFFDInfo(CGeometry *geometry, CConfig *config, CFreeFormDefBox **FFDBox, string val_mesh_filename, bool val_fullmesh);
	
	/*! 
	 * \brief Write the Free Form information in the SU2 file.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] FFDBox - Array with all the free forms FFDBoxs of the computation.
	 * \param[in] val_mesh_filename - Name of the grid output file.
	 */		
	void WriteFFDInfo(CGeometry *geometry, CConfig *config, string val_mesh_filename);
  
  /*!
	 * \brief Write the Free Form information in the SU2 file.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] FFDBox - Array with all the free forms FFDBoxs of the computation.
	 * \param[in] val_mesh_filename - Name of the grid output file.
	 */
	void WriteFFDInfo(CGeometry *geometry, CConfig *config, CFreeFormDefBox **FFDBox, string val_mesh_filename);
	
	/*! 
	 * \brief Get information about if there is a complete FFDBox definition, or it is necessary to 
	 *        compute the parametric coordinates.
	 * \return <code>TRUE</code> if the input grid file has a complete information; otherwise <code>FALSE</code>.
	 */		
	bool GetFFDBoxDefinition(void);
	
	/*! 
	 * \brief Obtain the number of FFDBoxs.
	 * \return Number of FFD FFDBoxs.
	 */		
	unsigned short GetnFFDBox(void);
	
	/*! 
	 * \brief Obtain the number of levels.
	 * \return Number of FFD levels.
	 */		
	unsigned short GetnLevel(void);
	
};

#include "grid_movement_structure.inl"
