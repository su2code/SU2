/*!
 * \file CFreeFormDefBox.hpp
 * \brief Headers of the CFreeFormDefBox class.
 * \author F. Palacios & A. Galdran.
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "CGridMovement.hpp"
#include "CFreeFormBlending.hpp"

/*!
 * \class CFreeFormDefBox
 * \brief Class for defining the free form FFDBox structure.
 * \author F. Palacios & A. Galdran.
 */
class CFreeFormDefBox : public CGridMovement {
 public:
  unsigned short nDim;                     /*!< \brief Number of dimensions of the problem. */
  unsigned short nCornerPoints,            /*!< \brief Number of corner points of the FFDBox. */
      nControlPoints, nControlPoints_Copy; /*!< \brief Number of control points of the FFDBox. */
  su2double **Coord_Corner_Points,         /*!< \brief Coordinates of the corner points. */
      ****Coord_Control_Points,            /*!< \brief Coordinates of the control points. */
      ****ParCoord_Control_Points,         /*!< \brief Coordinates of the control points. */
      ****Coord_Control_Points_Copy,       /*!< \brief Coordinates of the control points (copy). */
      ****Coord_SupportCP{nullptr};        /*!< \brief Coordinates of the support control points. */
  unsigned short lOrder, lOrder_Copy,      /*!< \brief Order of the FFDBox in the i direction. */
      mOrder, mOrder_Copy,                 /*!< \brief Order of the FFDBox in the j direction. */
      nOrder, nOrder_Copy;                 /*!< \brief Order of the FFDBox in the k direction. */
  unsigned short lDegree, lDegree_Copy,    /*!< \brief Degree of the FFDBox in the i direction. (lOrder - 1)*/
      mDegree, mDegree_Copy,               /*!< \brief Degree of the FFDBox in the j direction. (mOrder - 1)*/
      nDegree, nDegree_Copy;               /*!< \brief Degree of the FFDBox in the k direction. (nOrder - 1)*/
  su2double *ParamCoord, *ParamCoord_,     /*!< \brief Parametric coordinates of a point. */
      *cart_coord, *cart_coord_;           /*!< \brief Cartesian coordinates of a point. */
  su2double ObjFunc;                       /*!< \brief Objective function of the point inversion process. */
  su2double* Gradient;                     /*!< \brief Gradient of the point inversion process. */
  su2double** Hessian;                     /*!< \brief Hessian of the point inversion process. */
  su2double MaxCoord[3];                   /*!< \brief Maximum coordinates of the FFDBox. */
  su2double MinCoord[3];                   /*!< \brief Minimum coordinates of the FFDBox. */
  string Tag;                              /*!< \brief Tag to identify the FFDBox. */
  unsigned short Level;                    /*!< \brief Nested level of the FFD box. */

  vector<su2double> CartesianCoord[3];  /*!< \brief Vector with all the cartesian coordinates in the FFD FFDBox. */
  vector<su2double> ParametricCoord[3]; /*!< \brief Vector with all the parametrics coordinates in the FFD FFDBox. */
  vector<unsigned short> MarkerIndex;   /*!< \brief Vector with all markers in the FFD FFDBox. */
  vector<unsigned long> VertexIndex;    /*!< \brief Vector with all vertex index in the FFD FFDBox. */
  vector<unsigned long> PointIndex;     /*!< \brief Vector with all points index in the FFD FFDBox. */
  unsigned long nSurfacePoint;          /*!< \brief Number of surfaces in the FFD FFDBox. */
  vector<string> ParentFFDBox;          /*!< \brief Vector with all the parent FFD FFDBox. */
  vector<string> ChildFFDBox;           /*!< \brief Vector with all the child FFD FFDBox. */
  vector<unsigned short> Fix_IPlane;    /*!< \brief Fix FFD I plane. */
  vector<unsigned short> Fix_JPlane;    /*!< \brief Fix FFD J plane. */
  vector<unsigned short> Fix_KPlane;    /*!< \brief Fix FFD K plane. */

  CFreeFormBlending** BlendingFunction;

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
  CFreeFormDefBox(const unsigned short Degree[], unsigned short BSplineOrder[], unsigned short kind_blending);

  /*!
   * \brief Destructor of the class.
   */
  ~CFreeFormDefBox(void) override;

  /*!
   * \brief Define the I planes to to fix in a FFD box.
   * \param[in] val_plane - Index of the plane to fix.
   */
  inline void Set_Fix_IPlane(unsigned short val_plane) { Fix_IPlane.push_back(val_plane); }

  /*!
   * \brief Define the I planes to to fix in a FFD box.
   * \param[in] val_plane - Index of the plane to fix.
   */
  inline void Set_Fix_JPlane(unsigned short val_plane) { Fix_JPlane.push_back(val_plane); }

  /*!
   * \brief Define the I planes to to fix in a FFD box.
   * \param[in] val_plane - Index of the plane to fix.
   */
  inline void Set_Fix_KPlane(unsigned short val_plane) { Fix_KPlane.push_back(val_plane); }

  /*!
   * \brief Define the I planes to to fix in a FFD box.
   * \param[in] val_plane - Index of the plane to fix.
   */
  inline unsigned short Get_Fix_IPlane(unsigned short val_index) { return Fix_IPlane[val_index]; }

  /*!
   * \brief Define the I planes to to fix in a FFD box.
   * \param[in] val_plane - Index of the plane to fix.
   */
  inline unsigned short Get_Fix_JPlane(unsigned short val_index) { return Fix_JPlane[val_index]; }

  /*!
   * \brief Define the I planes to to fix in a FFD box.
   * \param[in] val_plane - Index of the plane to fix.
   */
  inline unsigned short Get_Fix_KPlane(unsigned short val_index) { return Fix_KPlane[val_index]; }

  /*!
   * \brief Define the I planes to to fix in a FFD box.
   * \param[in] val_plane - Index of the plane to fix.
   */
  inline unsigned short Get_nFix_IPlane(void) const { return Fix_IPlane.size(); }

  /*!
   * \brief Define the I planes to to fix in a FFD box.
   * \param[in] val_plane - Index of the plane to fix.
   */
  inline unsigned short Get_nFix_JPlane(void) const { return Fix_JPlane.size(); }

  /*!
   * \brief Define the I planes to to fix in a FFD box.
   * \param[in] val_plane - Index of the plane to fix.
   */
  inline unsigned short Get_nFix_KPlane(void) const { return Fix_KPlane.size(); }

  /*!
   * \brief Add to the vector of markers a new marker.
   * \param[in] val_iMarker - New marker inside the FFD box.
   */
  inline void Set_MarkerIndex(unsigned short val_iMarker) { MarkerIndex.push_back(val_iMarker); }

  /*!
   * \brief Add to the vector of vertices a new vertex.
   * \param[in] val_iVertex - New vertex inside the FFD box.
   */
  inline void Set_VertexIndex(unsigned long val_iVertex) { VertexIndex.push_back(val_iVertex); }

  /*!
   * \brief Add to the vector of points a new point.
   * \param[in] val_iPoint - New point inside the FFD box.
   */
  inline void Set_PointIndex(unsigned long val_iPoint) { PointIndex.push_back(val_iPoint); }

  /*!
   * \brief Add to the vector of cartesian coordinates a new coordinate.
   * \param[in] val_coord - New coordinate inside the FFD box.
   */
  inline void Set_CartesianCoord(su2double* val_coord) {
    CartesianCoord[0].push_back(val_coord[0]);
    CartesianCoord[1].push_back(val_coord[1]);
    CartesianCoord[2].push_back(val_coord[2]);
  }

  /*!
   * \brief Adds to the vector of cartesian coordinates.
   * \param[in] val_coord - New coord inside FFD box.
   * \param[in] val_iSurfacePoints - Surface points of FFD box.
   */
  inline void Set_CartesianCoord(const su2double* val_coord, unsigned long val_iSurfacePoints) {
    CartesianCoord[0][val_iSurfacePoints] = val_coord[0];
    CartesianCoord[1][val_iSurfacePoints] = val_coord[1];
    CartesianCoord[2][val_iSurfacePoints] = val_coord[2];
  }

  /*!
   * \brief Add to the vector of parametric coordinates a new coordinate.
   * \param[in] val_coord - New coordinate inside the FFD box.
   */
  inline void Set_ParametricCoord(su2double* val_coord) {
    ParametricCoord[0].push_back(val_coord[0]);
    ParametricCoord[1].push_back(val_coord[1]);
    ParametricCoord[2].push_back(val_coord[2]);
  }

  /*!
   * \brief Add to the vector of parent FFDBoxes a new FFD FFDBox.
   * \param[in] val_iParentFFDBox - New parent FFDBox in the vector.
   */
  inline void SetParentFFDBox(string val_iParentFFDBox) { ParentFFDBox.push_back(val_iParentFFDBox); }

  /*!
   * \brief Add to the vector of child FFDBoxes a new FFD FFDBox.
   * \param[in] val_iChildFFDBox - New child FFDBox in the vector.
   */
  inline void SetChildFFDBox(string val_iChildFFDBox) { ChildFFDBox.push_back(val_iChildFFDBox); }

  /*!
   * \brief Adds to the set of Parametric coordinates.
   * \param[in] val_coord - New coord inside FFD box.
   * \param[in] val_iSurfacePoints - Surface points of FFD box.
   */
  inline void Set_ParametricCoord(const su2double* val_coord, unsigned long val_iSurfacePoints) {
    ParametricCoord[0][val_iSurfacePoints] = val_coord[0];
    ParametricCoord[1][val_iSurfacePoints] = val_coord[1];
    ParametricCoord[2][val_iSurfacePoints] = val_coord[2];
  }

  /*!
   * \brief Get index of the marker.
   * \param[in] val_iSurfacePoints - Surface points of FFD box.
   */
  inline unsigned short Get_MarkerIndex(unsigned long val_iSurfacePoints) { return MarkerIndex[val_iSurfacePoints]; }

  /*!
   * \brief Get index of the marker.
   * \param[in] Get_VertexIndex - Surface points of FFD box.
   */
  inline unsigned long Get_VertexIndex(unsigned long val_iSurfacePoints) { return VertexIndex[val_iSurfacePoints]; }

  /*!
   * \brief Get index of the point.
   * \param[in] Get_VertexIndex - Surface points of FFD box.
   */
  inline unsigned long Get_PointIndex(unsigned long val_iSurfacePoints) { return PointIndex[val_iSurfacePoints]; }

  /*!
   * \brief Get Cartesian coordinates.
   * \param[in] Get_VertexIndex - Surface points of FFD box.
   */
  inline su2double* Get_CartesianCoord(unsigned long val_iSurfacePoints) {
    cart_coord_[0] = CartesianCoord[0][val_iSurfacePoints];
    cart_coord_[1] = CartesianCoord[1][val_iSurfacePoints];
    cart_coord_[2] = CartesianCoord[2][val_iSurfacePoints];
    return cart_coord_;
  }

  /*!
   * \brief Get parametric coordinates.
   * \param[in] Get_VertexIndex - Surface points of FFD box.
   */
  inline su2double* Get_ParametricCoord(unsigned long val_iSurfacePoints) {
    ParamCoord_[0] = ParametricCoord[0][val_iSurfacePoints];
    ParamCoord_[1] = ParametricCoord[1][val_iSurfacePoints];
    ParamCoord_[2] = ParametricCoord[2][val_iSurfacePoints];
    return ParamCoord_;
  }

  /*!
   * \brief Get number of surface points.
   */
  inline unsigned long GetnSurfacePoint(void) const { return PointIndex.size(); }

  /*!
   * \brief Get number of parent FFD boxes.
   */
  inline unsigned short GetnParentFFDBox(void) const { return ParentFFDBox.size(); }

  /*!
   * \brief Get number of child FFD boxes.
   */
  inline unsigned short GetnChildFFDBox(void) const { return ChildFFDBox.size(); }

  /*!
   * \brief Get tag of parent FFD box.
   * \param[in] val_ParentFFDBox - idex of parent FFD box.
   */
  inline string GetParentFFDBoxTag(unsigned short val_ParentFFDBox) { return ParentFFDBox[val_ParentFFDBox]; }

  /*!
   * \brief Get tag of child FFD box.
   * \param[in] val_ChildFFDBox - index of child FFD box.
   */
  inline string GetChildFFDBoxTag(unsigned short val_ChildFFDBox) { return ChildFFDBox[val_ChildFFDBox]; }

  /*!
   * \brief Change the the position of the corners of the unitary FFDBox,
   *        and find the position of the control points for the FFDBox
   * \param[in] FFDBox - Original FFDBox where we want to compute the control points.
   */
  void SetSupportCPChange(CFreeFormDefBox* FFDBox);

  /*!
   * \brief Set the number of corner points.
   * \param[in] val_ncornerpoints - Number of corner points.
   */
  inline void SetnCornerPoints(unsigned short val_ncornerpoints) { nCornerPoints = val_ncornerpoints; }

  /*!
   * \brief Get the number of corner points.
   * \return Number of corner points.
   */
  inline unsigned short GetnCornerPoints(void) const { return nCornerPoints; }

  /*!
   * \brief Get the number of control points.
   * \return Number of control points.
   */
  inline unsigned short GetnControlPoints(void) const { return nControlPoints; }

  /*!
   * \brief Get the number of control points.
   * \return Number of control points.
   */
  inline void SetnControlPoints(void) { nControlPoints = lOrder * mOrder * nOrder; }

  /*!
   * \brief Get the number of numerical points on the surface.
   * \return Number of numerical points on the surface.
   */
  inline unsigned long GetnSurfacePoints(void) { return 0; }

  /*!
   * \brief Set the corner point for the unitary FFDBox.
   */
  void SetUnitCornerPoints(void);

  /*!
   * \brief Set the coordinates of the corner points.
   * \param[in] val_coord - Coordinates of the corner point with index <i>val_icornerpoints</i>.
   * \param[in] val_icornerpoints - Index of the corner point.
   */
  inline void SetCoordCornerPoints(const su2double* val_coord, unsigned short val_icornerpoints) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) Coord_Corner_Points[val_icornerpoints][iDim] = val_coord[iDim];
  }

  /*!
   * \overload
   * \param[in] val_xcoord - X coordinate of the corner point with index <i>val_icornerpoints</i>.
   * \param[in] val_ycoord - Y coordinate of the corner point with index <i>val_icornerpoints</i>.
   * \param[in] val_zcoord - Z coordinate of the corner point with index <i>val_icornerpoints</i>.
   * \param[in] val_icornerpoints - Index of the corner point.
   */
  inline void SetCoordCornerPoints(su2double val_xcoord, su2double val_ycoord, su2double val_zcoord,
                                   unsigned short val_icornerpoints) {
    Coord_Corner_Points[val_icornerpoints][0] = val_xcoord;
    Coord_Corner_Points[val_icornerpoints][1] = val_ycoord;
    Coord_Corner_Points[val_icornerpoints][2] = val_zcoord;
  }

  /*!
   * \brief Set the coordinates of the control points.
   * \param[in] val_coord - Coordinates of the control point.
   * \param[in] iDegree - Index of the FFDBox, i direction.
   * \param[in] jDegree - Index of the FFDBox, j direction.
   * \param[in] kDegree - Index of the FFDBox, k direction.
   */
  inline void SetCoordControlPoints(const su2double* val_coord, unsigned short iDegree, unsigned short jDegree,
                                    unsigned short kDegree) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Coord_Control_Points[iDegree][jDegree][kDegree][iDim] = val_coord[iDim];
    }
  }

  /*!
   * \brief Set the coordinates of the control points.
   * \param[in] val_coord - Coordinates of the control point.
   * \param[in] iDegree - Index of the FFDBox, i direction.
   * \param[in] jDegree - Index of the FFDBox, j direction.
   * \param[in] kDegree - Index of the FFDBox, k direction.
   */
  inline void SetCoordControlPoints_Copy(const su2double* val_coord, unsigned short iDegree, unsigned short jDegree,
                                         unsigned short kDegree) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Coord_Control_Points_Copy[iDegree][jDegree][kDegree][iDim] = val_coord[iDim];
    }
  }

  /*!
   * \brief Set the coordinates of the control points.
   * \param[in] val_coord - Coordinates of the control point.
   * \param[in] iDegree - Index of the FFDBox, i direction.
   * \param[in] jDegree - Index of the FFDBox, j direction.
   * \param[in] kDegree - Index of the FFDBox, k direction.
   */
  inline void SetParCoordControlPoints(const su2double* val_coord, unsigned short iDegree, unsigned short jDegree,
                                       unsigned short kDegree) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      ParCoord_Control_Points[iDegree][jDegree][kDegree][iDim] = val_coord[iDim];
  }

  /*!
   * \brief Get the coordinates of the corner points.
   * \param[in] val_dim - Index of the coordinate (x, y, z).
   * \param[in] val_icornerpoints - Index of the corner point.
   * \return Coordinate <i>val_dim</i> of the corner point <i>val_icornerpoints</i>.
   */
  inline su2double GetCoordCornerPoints(unsigned short val_dim, unsigned short val_icornerpoints) const {
    return Coord_Corner_Points[val_icornerpoints][val_dim];
  }

  /*!
   * \brief Get the coordinates of the corner points.
   * \param[in] val_icornerpoints - Index of the corner point.
   * \return Pointer to the coordinate vector of the corner point <i>val_icornerpoints</i>.
   */
  inline su2double* GetCoordCornerPoints(unsigned short val_icornerpoints) const {
    return Coord_Corner_Points[val_icornerpoints];
  }

  /*!
   * \brief Get the coordinates of the control point.
   * \param[in] val_iindex - Value of the local i index of the control point.
   * \param[in] val_jindex - Value of the local j index of the control point.
   * \param[in] val_kindex - Value of the local k index of the control point.
   * \return Pointer to the coordinate vector of the control point with local index (i, j, k).
   */
  inline su2double* GetCoordControlPoints(unsigned short val_iindex, unsigned short val_jindex,
                                          unsigned short val_kindex) const {
    return Coord_Control_Points[val_iindex][val_jindex][val_kindex];
  }

  /*!
   * \brief Get the parametric coordinates of the control point.
   * \param[in] val_iindex - Value of the local i index of the control point.
   * \param[in] val_jindex - Value of the local j index of the control point.
   * \param[in] val_kindex - Value of the local k index of the control point.
   * \return Pointer to the coordinate vector of the control point with local index (i, j, k).
   */
  inline su2double* GetParCoordControlPoints(unsigned short val_iindex, unsigned short val_jindex,
                                             unsigned short val_kindex) const {
    return ParCoord_Control_Points[val_iindex][val_jindex][val_kindex];
  }

  /*!
   * \brief Set the control points in a parallelepiped (hexahedron).
   */
  void SetControlPoints_Parallelepiped(void);

  /*!
   * \brief Set the control points of the final chuck in a unitary hexahedron free form.
   * \param[in] FFDBox - Original FFDBox where we want to compute the control points.
   */
  void SetSupportCP(CFreeFormDefBox* FFDBox);

  /*!
   * \brief Set the new value of the coordinates of the control points.
   * \param[in] val_index - Local index (i, j, k) of the control point.
   * \param[in] movement - Movement of the control point.
   */
  inline void SetControlPoints(const unsigned short* val_index, const su2double* movement) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Coord_Control_Points[val_index[0]][val_index[1]][val_index[2]][iDim] += movement[iDim];
  }

  /*!
   * \brief Set the original value of the control points.
   */
  inline void SetOriginalControlPoints() {
    for (unsigned short iDegree = 0; iDegree <= lDegree_Copy; iDegree++)
      for (unsigned short jDegree = 0; jDegree <= mDegree_Copy; jDegree++)
        for (unsigned short kDegree = 0; kDegree <= nDegree_Copy; kDegree++)
          for (unsigned short iDim = 0; iDim < nDim; iDim++)
            Coord_Control_Points[iDegree][jDegree][kDegree][iDim] =
                Coord_Control_Points_Copy[iDegree][jDegree][kDegree][iDim];

    lDegree = lDegree_Copy;
    mDegree = mDegree_Copy;
    nDegree = nDegree_Copy;
    lOrder = lOrder_Copy;
    mOrder = mOrder_Copy;
    nOrder = nOrder_Copy;
    nControlPoints = nControlPoints_Copy;
  }

  /*!
   * \brief Set the tecplot file of the FFD chuck structure.
   * \param[in] iFFDBox - Index of the FFD box.
   * \param[in] original - Original box (before deformation).
   */
  void SetTecplot(CGeometry* geometry, unsigned short iFFDBox, bool original);

  /*!
   * \brief Set the paraview file of the FFD chuck structure.
   * \param[in] iFFDBox - Index of the FFD box.
   * \param[in] original - Original box (before deformation).
   */
  void SetParaview(CGeometry* geometry, unsigned short iFFDBox, bool original);

  /*!
   * \brief Set the CGNS file of the FFD chuck structure.
   * \param[in] iFFDBox - Index of the FFD box.
   * \param[in] original - Original box (before deformation).
   */
  void SetCGNS(CGeometry* geometry, unsigned short iFFDBox, bool original);

  /*!
   * \brief Set Cylindrical to Cartesians_ControlPoints.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCyl2Cart_ControlPoints(CConfig* config);

  /*!
   * \brief Set Cartesians to Cylindrical ControlPoints.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCart2Cyl_ControlPoints(CConfig* config);

  /*!
   * \brief Set Cylindrical to Cartesians_CornerPoints.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCyl2Cart_CornerPoints(CConfig* config);

  /*!
   * \brief Set Cartesians to Cylindrical CornerPoints.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCart2Cyl_CornerPoints(CConfig* config);

  /*!
   * \brief Set Spherical to Cartesians ControlPoints.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSphe2Cart_ControlPoints(CConfig* config);

  /*!
   * \brief SetCartesians to Spherical ControlPoints.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCart2Sphe_ControlPoints(CConfig* config);

  /*!
   * \brief Set Spherical to Cartesians_CornerPoints.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSphe2Cart_CornerPoints(CConfig* config);

  /*!
   * \brief Set Cartesians to Spherical Corner Points.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCart2Sphe_CornerPoints(CConfig* config);

  /*!
   * \brief Set the cartesian coords of a point in R^3 and convert them to the parametric coords of
   *        our parametrization of a paralellepiped.
   * \param[in] cart_coord - Cartesian coordinates of a point.
   * \return Pointer to the parametric coordinates of a point.
   */
  su2double* GetParametricCoord_Analytical(const su2double* cart_coord);

  /*!
   * \brief Iterative strategy for computing the parametric coordinates.
   * \param[in] xyz - Cartesians coordinates of the target point.
   * \param[in] guess - Initial guess for doing the parametric coordinates search.
   * \param[in] tol - Level of convergence of the iterative method.
   * \param[in] it_max - Maximal number of iterations.
   * \return Parametric coordinates of the point.
   */
  su2double* GetParametricCoord_Iterative(unsigned long iPoint, su2double* xyz, const su2double* guess,
                                          CConfig* config);

  /*!
   * \brief Compute the cross product.
   * \param[in] v1 - First input vector.
   * \param[in] v2 - Second input vector.
   * \param[out] v3 - Output vector wuth the cross product.
   */
  inline void CrossProduct(const su2double* v1, const su2double* v2, su2double* v3) {
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }

  /*!
   * \brief Compute the doc product.
   * \param[in] v1 - First input vector.
   * \param[in] v2 - Sencond input vector.
   * \return Dot product between <i>v1</i>, and <i>v2</i>.
   */
  inline su2double DotProduct(const su2double* v1, const su2double* v2) {
    su2double scalar = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    return scalar;
  }

  /*!
   * \brief Here we take the parametric coords of a point in the box and we convert them to the
   *        physical cartesian coords by plugging the ParamCoords on the Bezier parameterization of our box.
   * \param[in] ParamCoord - Parametric coordinates of a point.
   * \return Pointer to the cartesian coordinates of a point.
   */
  su2double* EvalCartesianCoord(su2double* ParamCoord) const;

  /*!
   * \brief Get the order in the l direction of the FFD FFDBox.
   * \return Order in the l direction of the FFD FFDBox.
   */
  inline unsigned short GetlOrder(void) const { return lOrder; }

  /*!
   * \brief Get the order in the m direction of the FFD FFDBox.
   * \return Order in the m direction of the FFD FFDBox.
   */
  inline unsigned short GetmOrder(void) const { return mOrder; }

  /*!
   * \brief Get the order in the n direction of the FFD FFDBox.
   * \return Order in the n direction of the FFD FFDBox.
   */
  inline unsigned short GetnOrder(void) const { return nOrder; }

  /*!
   * \brief Get the order in the l direction of the FFD FFDBox.
   * \return Order in the l direction of the FFD FFDBox.
   */
  inline void SetlOrder(unsigned short val_lOrder) {
    lOrder = val_lOrder;
    lDegree = lOrder - 1;
  }

  /*!
   * \brief Get the order in the m direction of the FFD FFDBox.
   * \return Order in the m direction of the FFD FFDBox.
   */
  inline void SetmOrder(unsigned short val_mOrder) {
    mOrder = val_mOrder;
    mDegree = mOrder - 1;
  }

  /*!
   * \brief Get the order in the n direction of the FFD FFDBox.
   * \return Order in the n direction of the FFD FFDBox.
   */
  inline void SetnOrder(unsigned short val_nOrder) {
    nOrder = val_nOrder;
    nDegree = nOrder - 1;
  }

  /*!
   * \brief Returns true if the point is inside the FFD.
   * \param[in] coord - Coordinate of the point to check.
   */
  bool CheckPointInsideFFD(const su2double* coord) const;

  /*!
   * \brief Set the zone of the computational domain that is going to be deformed.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iFFDBox - Index of the FFDBox.
   */
  // this routine is not used. We should consider deleting it.
  void SetDeformationZone(CGeometry* geometry, CConfig* config, unsigned short iFFDBox) const;

  /*!
   * \brief The routine computes the gradient of F(u, v, w) = ||X(u, v, w)-(x, y, z)||^2  evaluated at (u, v, w).
   * \param[in] val_coord - Parametric coordiates of the target point.
   * \param[in] xyz - Cartesians coordinates of the point.
   * \param[in] analytical - Compute the analytical gradient.
   * \return Value of the analytical gradient.
   */
  su2double* GetFFDGradient(su2double* val_coord, su2double* xyz);

  /*!
   * \brief The routine that computes the Hessian of F(u, v, w) = ||X(u, v, w)-(x, y, z)||^2 evaluated at (u, v, w)
   *        Input: (u, v, w), (x, y, z)
   *        Output: Hessian F (u, v, w).
   * \param[in] uvw - Current value of the parametrics coordinates.
   * \param[in] xyz - Cartesians coordinates of the target point to compose the functional.
   * \param[in] val_Hessian - Value of the hessian.
   */
  void GetFFDHessian(su2double* uvw, su2double* xyz, su2double** val_Hessian);

  /*!
   * \brief An auxiliary routine to help us compute the gradient of F(u, v, w) = ||X(u, v, w)-(x, y, z)||^2 =
   *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)^2+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z)^2
   *        Input: val_t, val_diff (to identify the index of the Bernstein polynomail we differentiate), the i, j, k ,
   * l, m, n E.G.: val_diff=2 => we differentiate w.r.t. w  (val_diff=0,1, or 2) Output: d [B_i^l*B_j^m *B_k^n] / d
   * val_diff (val_u, val_v, val_w). \param[in] uvw - __________. \param[in] val_diff - __________. \param[in] ijk -
   * __________. \param[in] lmn - Degree of the FFD box. \return __________.
   */
  su2double GetDerivative1(su2double* uvw, unsigned short val_diff, unsigned short* ijk, unsigned short* lmn) const;

  /*!
   * \brief An auxiliary routine to help us compute the gradient of F(u, v, w) = ||X(u, v, w)-(x, y, z)||^2 =
   *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)^2+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z)^2
   *        Input: (u, v, w), dim , xyz=(x, y, z), l, m, n E.G.: dim=2 => we use the third coordinate of the control
   * points, and the z-coordinate of xyz  (0<=dim<=2) Output: 2* ( (Sum_{i, j, k}^l, m, n P_{ijk}[dim] B_i^l[u] B_j^m[v]
   * B_k^n[w]) - xyz[dim]). \param[in] uvw - __________. \param[in] dim - __________. \param[in] xyz - __________.
   * \param[in] lmn - Degree of the FFD box.
   * \return __________.
   */
  su2double GetDerivative2(su2double* uvw, unsigned short dim, const su2double* xyz, const unsigned short* lmn) const;

  /*!
   * \brief An auxiliary routine to help us compute the gradient of F(u, v, w) = ||X(u, v, w)-(x, y, z)||^2 =
   *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z)
   * \param[in] uvw - Parametric coordiates of the point.
   * \param[in] dim - Value of the coordinate to be differentiate.
   * \param[in] diff_this - Diferentiation with respect this coordinate.
   * \param[in] lmn - Degree of the FFD box.
   * \return Sum_{i, j, k}^{l, m, n} [one of them with -1,
   *        depending on diff_this=0,1 or 2] P_{ijk}[dim] * (B_i^l[u] B_j^m[v] B_k^n[w])--one of them diffrentiated;
   *        which? diff_thiss will tell us ; E.G.: dim=2, diff_this=1 => we use the third coordinate of the control
   *        points, and derivate de v-Bersntein polynomial (use m-1 when summing!!).
   */
  su2double GetDerivative3(su2double* uvw, unsigned short dim, unsigned short diff_this, unsigned short* lmn) const;

  /*!
   * \brief An auxiliary routine to help us compute the Hessian of F(u, v, w) = ||X(u, v, w)-(x, y, z)||^2 =
   *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z)
   *        Input: val_t, val_diff, val_diff2 (to identify the index of the Bernstein polynomials we differentiate), the
   * i, j, k , l, m, n E.G.: val_diff=1, val_diff2=2  =>  we differentiate w.r.t. v and w  (val_diff=0,1, or 2) E.G.:
   * val_diff=0, val_diff2=0 => we differentiate w.r.t. u two times Output: [d [B_i^l*B_j^m *B_k^n]/d val_diff *d
   * [B_i^l*B_j^m *B_k^n]/d val_diff2] (val_u, val_v, val_w) . \param[in] uvw - __________. \param[in] val_diff -
   * __________. \param[in] val_diff2 - __________. \param[in] ijk - __________. \param[in] lmn - Degree of the FFD box.
   * \return __________.
   */
  su2double GetDerivative4(su2double* uvw, unsigned short val_diff, unsigned short val_diff2, unsigned short* ijk,
                           unsigned short* lmn) const;

  /*!
   * \brief An auxiliary routine to help us compute the Hessian of F(u, v, w) = ||X(u, v, w)-(x, y, z)||^2 =
   *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z)
   *        Input: (u, v, w), dim , diff_this, diff_this_also, xyz=(x, y, z), l, m, n
   *        Output:
   *        Sum_{i, j, k}^{l, m, n} [two of them with -1, depending on diff_this, diff_this_also=0,1 or 2]
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
  su2double GetDerivative5(su2double* uvw, unsigned short dim, unsigned short diff_this, unsigned short diff_this_also,
                           unsigned short* lmn) const;

  /*!
   * \brief Euclidean norm of a vector.
   * \param[in] a - _______.
   * \return __________.
   */
  inline su2double GetNorm(const su2double* a) {
    su2double norm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    return norm;
  }

  /*!
   * \brief Set the tag that identify a FFDBox.
   * \param[in] val_tag - value of the tag.
   */
  inline void SetTag(string val_tag) { Tag = val_tag; }

  /*!
   * \brief Get the tag that identify a FFDBox.
   * \return Value of the tag that identigy the FFDBox.
   */
  inline string GetTag() const { return Tag; }

  /*!
   * \brief Set the nested level of the FFDBox.
   * \param[in] val_level - value of the level.
   */
  inline void SetLevel(unsigned short val_level) { Level = val_level; }

  /*!
   * \brief Get the nested level of the FFDBox.
   * \return Value of the nested level of the the FFDBox.
   */
  inline unsigned short GetLevel() const { return Level; }
};
