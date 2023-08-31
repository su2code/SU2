/*!
 * \file CADTElemClass.hpp
 * \brief Class for storing an ADT of (linear) elements in an arbitrary number of dimensions.
 * \author E. van der Weide
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

#include "./CADTBaseClass.hpp"
#include "./CBBoxTargetClass.hpp"
#include "../parallelization/omp_structure.hpp"

/*!
 * \class CADTElemClass
 * \ingroup ADT
 * \brief  Class for storing an ADT of (linear) elements in an arbitrary number of dimensions.
 * \author E. van der Weide
 * \version 8.0.0 "Harrier"
 */
class CADTElemClass : public CADTBaseClass {
 private:
  unsigned short nDim; /*!< \brief Number of spatial dimensions. */

  vector<su2double> coorPoints; /*!< \brief Vector, which contains the coordinates
                                            of the points in the ADT. */
  vector<su2double> BBoxCoor;   /*!< \brief Vector, which contains the coordinates
                                            of the bounding boxes of the elements. */

  vector<unsigned short> elemVTK_Type; /*!< \brief Vector, which the type of the elements
                                                   using the VTK convention. */
  vector<unsigned long> nDOFsPerElem;  /*!< \brief Vector, which contains the number of DOFs
                                                   of the elements in the ADT in cumulative
                                                   storage format. */
  vector<unsigned long> elemConns;     /*!< \brief Vector, which contains the connectivities
                                                   of the elements in the ADT. */
  vector<unsigned short> localMarkers; /*!< \brief Vector, which contains the marker ID's
                                                   of the elements in the ADT. */
  vector<unsigned long> localElemIDs;  /*!< \brief Vector, which contains the local element ID's
                                                   of the elements in the ADT. */
  vector<int> ranksOfElems;            /*!< \brief Vector, which contains the ranks
                                                   of the elements in the ADT. */
#ifdef HAVE_OMP
  vector<vector<CBBoxTargetClass> > BBoxTargets; /*!< \brief Vector, used to store possible bounding box
                                                             candidates during the nearest element search. */
#else
  array<vector<CBBoxTargetClass>, 1> BBoxTargets;
#endif
 public:
  /*!
   * \brief Constructor of the class.
   * \param[in]     val_nDim     Number of spatial dimensions of the problem.
   * \param[in]     val_coor     Coordinates of the local points to be stored in the ADT.
   * \param[in,out] val_connElem Local connectivity of the elements to be stored in the ADT.
                                 In parallel mode the connectivities are corrected for the
                                 offset in node numbers of the rank.
   * \param[in]     val_VTKElem  Type of the elements using the VTK convention.
   * \param[in]     val_markerID Markers of the local elements.
   * \param[in]     val_elemID   Local element IDs of the elements.
   * \param[in]     globalTree   Whether or not a global tree must be built. If false
                                 a local ADT is built.
   */
  CADTElemClass(unsigned short val_nDim, vector<su2double>& val_coor, vector<unsigned long>& val_connElem,
                vector<unsigned short>& val_VTKElem, vector<unsigned short>& val_markerID,
                vector<unsigned long>& val_elemID, const bool globalTree);

  /*!
   * \brief Function, which determines the element that contains the given coordinate.
   * \note This simply forwards the call to the implementation function selecting the right
   *       working variables for the current thread.
   * \param[in]  coor             Coordinate which the element must contain.
   * \param[out] markerID         Local marker ID of the element containing the coordinate.
   * \param[out] elemID           Local element ID of the element containing the coordinate.
   * \param[out] rankID           Rank on which element containing the coordinate is stored.
   * \param[out] parCoor          Parametric coordinates of coor inside the element,
                                  which contains the coordinate.
   * \param[out] weightsInterpol  Interpolation weights of coor inside the element,
                                  which contains the coordinate.
   * \return                      True if an element is found, false if not.
   */
  inline bool DetermineContainingElement(const su2double* coor, unsigned short& markerID, unsigned long& elemID,
                                         int& rankID, su2double* parCoor, su2double* weightsInterpol) {
    const auto iThread = omp_get_thread_num();
    return DetermineContainingElement_impl(FrontLeaves[iThread], FrontLeavesNew[iThread], coor, markerID, elemID,
                                           rankID, parCoor, weightsInterpol);
  }

  /*!
   * \brief Function, which determines the nearest element in the ADT for the given coordinate.
   * \note This simply forwards the call to the implementation function selecting the right
   *       working variables for the current thread.
   * \param[in]  coor     Coordinate for which the nearest element in the ADT must be determined.
   * \param[out] dist     Distance to the nearest element in the ADT.
   * \param[out] markerID Local marker ID of the nearest element in the ADT.
   * \param[out] elemID   Local element ID of the nearest element in the ADT.
   * \param[out] rankID   Rank on which the nearest element in the ADT is stored.
   */
  inline void DetermineNearestElement(const su2double* coor, su2double& dist, unsigned short& markerID,
                                      unsigned long& elemID, int& rankID) {
    const auto iThread = omp_get_thread_num();
    DetermineNearestElement_impl(BBoxTargets[iThread], FrontLeaves[iThread], FrontLeavesNew[iThread], coor, dist,
                                 markerID, elemID, rankID);
  }

 private:
  /*!
   * \brief Implementation of DetermineContainingElement.
   * \note Working variables (first two) passed explicitly for thread safety.
   */
  bool DetermineContainingElement_impl(vector<unsigned long>& frontLeaves, vector<unsigned long>& frontLeavesNew,
                                       const su2double* coor, unsigned short& markerID, unsigned long& elemID,
                                       int& rankID, su2double* parCoor, su2double* weightsInterpol) const;

  /*!
   * \brief Implementation of DetermineNearestElement.
   * \note Working variables (first three) passed explicitly for thread safety.
   */
  void DetermineNearestElement_impl(vector<CBBoxTargetClass>& BBoxTargets, vector<unsigned long>& frontLeaves,
                                    vector<unsigned long>& frontLeavesNew, const su2double* coor, su2double& dist,
                                    unsigned short& markerID, unsigned long& elemID, int& rankID) const;

  /*!
   * \brief Function, which checks whether or not the given coordinate is
            inside the given element.
   * \param[in]  elemID          ID of the element for which the containment
                                 must be checked.
   * \param[in]  coor            Coordinate for which the containment must be checked.
   * \param[out] parCoor         Parametric coordinates of coor if it is inside the
                                 given element.
   * \param[out] weightsInterpol Interpolation weights of coor if it is inside the
                                 given element.
   * \return                     True if coor is inside the element and false otherwise.
   */
  bool CoorInElement(const unsigned long elemID, const su2double* coor, su2double* parCoor,
                     su2double* weightsInterpol) const;

  /*!
   * \brief Function, which checks whether or not the given coordinate is
            inside the given quadrilateral.
   * \param[in]  elemID          ID of the quadrilateral for which the containment
                                 must be checked.
   * \param[in]  coor            Coordinate for which the containment must be checked.
   * \param[out] parCoor         Parametric coordinates of coor if it is inside the
                                 given quadrilateral.
   * \param[out] weightsInterpol Interpolation weights of coor if it is inside the
                                 given quadrilateral.
   * \return                     True if coor is inside the quadrilateral and false otherwise.
   */
  bool CoorInQuadrilateral(const unsigned long elemID, const su2double* coor, su2double* parCoor,
                           su2double* weightsInterpol) const;

  /*!
   * \brief Function, which checks whether or not the given coordinate is
            inside the given triangle.
   * \param[in]  elemID          ID of the triangle for which the containment
                                 must be checked.
   * \param[in]  coor            Coordinate for which the containment must be checked.
   * \param[out] parCoor         Parametric coordinates of coor if it is inside the
                                 given triangle.
   * \param[out] weightsInterpol Interpolation weights of coor if it is inside the
                                 given triangle.
   * \return                     True if coor is inside the triangle and false otherwise.
   */
  bool CoorInTriangle(const unsigned long elemID, const su2double* coor, su2double* parCoor,
                      su2double* weightsInterpol) const;

  /*!
   * \brief Function, which checks whether or not the given coordinate is
            inside the given hexahedron.
   * \param[in]  elemID          ID of the hexahedron for which the containment
                                 must be checked.
   * \param[in]  coor            Coordinate for which the containment must be checked.
   * \param[out] parCoor         Parametric coordinates of coor if it is inside the
                                 given hexahedron.
   * \param[out] weightsInterpol Interpolation weights of coor if it is inside the
                                 given hexahedron.
   * \return                     True if coor is inside the hexahedron and false otherwise.
   */
  bool CoorInHexahedron(const unsigned long elemID, const su2double* coor, su2double* parCoor,
                        su2double* weightsInterpol) const;

  /*!
   * \brief Function, which checks whether or not the given coordinate is
            inside the given prism.
   * \param[in]  elemID          ID of the prism for which the containment
                                 must be checked.
   * \param[in]  coor            Coordinate for which the containment must be checked.
   * \param[out] parCoor         Parametric coordinates of coor if it is inside the
                                 given prism.
   * \param[out] weightsInterpol Interpolation weights of coor if it is inside the
                                 given prism.
   * \return                     True if coor is inside the prism and false otherwise.
   */
  bool CoorInPrism(const unsigned long elemID, const su2double* coor, su2double* parCoor,
                   su2double* weightsInterpol) const;

  /*!
   * \brief Function, which checks whether or not the given coordinate is
            inside the given pyramid.
   * \param[in]  elemID          ID of the pyramid for which the containment
                                 must be checked.
   * \param[in]  coor            Coordinate for which the containment must be checked.
   * \param[out] parCoor         Parametric coordinates of coor if it is inside the
                                 given pyramid.
   * \param[out] weightsInterpol Interpolation weights of coor if it is inside the
                                 given pyramid.
   * \return                     True if coor is inside the pyramid and false otherwise.
   */
  bool CoorInPyramid(const unsigned long elemID, const su2double* coor, su2double* parCoor,
                     su2double* weightsInterpol) const;

  /*!
   * \brief Function, which checks whether or not the given coordinate is
            inside the given tetrahedron.
   * \param[in]  elemID          ID of the tetrahedron for which the containment
                                 must be checked.
   * \param[in]  coor            Coordinate for which the containment must be checked.
   * \param[out] parCoor         Parametric coordinates of coor if it is inside the
                                 given tetrahedron.
   * \param[out] weightsInterpol Interpolation weights of coor if it is inside the
                                 given tetrahedron.
   * \return                     True if coor is inside the tetrahedron and false otherwise.
   */
  bool CoorInTetrahedron(const unsigned long elemID, const su2double* coor, su2double* parCoor,
                         su2double* weightsInterpol) const;

  /*!
   * \brief Function, which provides an initial guess for the parametric coordinates
            of the given point inside a hexahedron by splitting it into tetrahedra.
   * \param[in]  xRelC   Coordinates of the point to be investigated relative to
                         vertex 0 of the hexahedron.
   * \param[in]  xRel    Coordinates of the vertices of the hexahedron relative to
                         vertex 0.
   * \param[out] parCoor Initial guess of the parametric coordinates.
   * \return             True if the initial guess is within the hexahedron and
                         false otherwise.
   */
  bool InitialGuessContainmentHexahedron(const su2double xRelC[3], const su2double xRel[8][3],
                                         su2double* parCoor) const;

  /*!
   * \brief Function, which provides an initial guess for the parametric coordinates
            of the given point inside a prism by splitting it into tetrahedra.
   * \param[in]  xRelC   Coordinates of the point to be investigated relative to
                         vertex 0 of the prism.
   * \param[in]  xRel    Coordinates of the vertices of the prism relative to
                         vertex 0.
   * \param[out] parCoor Initial guess of the parametric coordinates.
   * \return             True if the initial guess is within the prism and
                         false otherwise.
   */
  bool InitialGuessContainmentPrism(const su2double xRelC[3], const su2double xRel[6][3], su2double* parCoor) const;

  /*!
   * \brief Function, which provides an initial guess for the parametric coordinates
            of the given point inside a pyramid by splitting it into tetrahedra.
   * \param[in]  xRelC   Coordinates of the point to be investigated relative to
                         vertex 0 of the pyramid.
   * \param[in]  xRel    Coordinates of the vertices of the pyramid relative to
                         vertex 0.
   * \param[out] parCoor Initial guess of the parametric coordinates.
   * \return             True if the initial guess is within the pyramid and
                         false otherwise.
   */
  bool InitialGuessContainmentPyramid(const su2double xRelC[3], const su2double xRel[5][3], su2double* parCoor) const;

  /*!
   * \brief Function, which computes the distance squared of the given coordinate
            to the given element.
   * \param[in]  elemID    ID of the element to which the distance must be determined.
   * \param[in]  coor      Coordinate for which the distance to the element must be determined.
   * \param[out] dist2Elem Distance squared from the coordinate to the element.
   */
  void Dist2ToElement(const unsigned long elemID, const su2double* coor, su2double& dist2Elem) const;
  /*!
   * \brief Function, which computes the distance squared of the given coordinate
            to a linear line element.
   * \param[in]  i0        Starting index in coorPoints, where the coordinates of the
                           first point of the line are stored.
   * \param[in]  i1        Starting index in coorPoints, where the coordinates of the
                           second point of the line are stored.
   * \param[in]  coor      Coordinate for which the distance to the line must be determined.
   * \param[out] dist2Line Distance squared from the coordinate to the line.
   */
  void Dist2ToLine(const unsigned long i0, const unsigned long i1, const su2double* coor, su2double& dist2Line) const;
  /*!
   * \brief Function, which computes the distance squared of the given coordinate
            to a linear quadrilateral element if the projection is inside the quad.
   * \param[in]  i0        Starting index in coorPoints, where the coordinates of the
                           first point of the quadrilateral are stored.
   * \param[in]  i1        Starting index in coorPoints, where the coordinates of the
                           second point of the quadrilateral are stored.
   * \param[in]  i2        Starting index in coorPoints, where the coordinates of the
                           third point of the quadrilateral are stored.
   * \param[in]  i3        Starting index in coorPoints, where the coordinates of the
                           fourth point of the quadrilateral are stored.
   * \param[in]  coor      Coordinate for which the distance to the quadrilateral
                           must be determined.
   * \param[out] r         Parametric coordinate of the projection.
   * \param[out] s         Parametric coordinate of the projection.
   * \param[out] dist2Quad Distance squared from the coordinate to the quadrilateral.
   * \return     True if the projection is inside the quadrilateral and false otherwise.
   */
  bool Dist2ToQuadrilateral(const unsigned long i0, const unsigned long i1, const unsigned long i2,
                            const unsigned long i3, const su2double* coor, su2double& r, su2double& s,
                            su2double& dist2Quad) const;
  /*!
   * \brief Function, which computes the distance squared of the given coordinate
            to a linear triangular element if the projection is inside the triangle.
   * \param[in]  i0        Starting index in coorPoints, where the coordinates of the
                           first point of the triangle are stored.
   * \param[in]  i1        Starting index in coorPoints, where the coordinates of the
                           second point of the triangle are stored.
   * \param[in]  i2        Starting index in coorPoints, where the coordinates of the
                           third point of the triangle are stored.
   * \param[in]  coor      Coordinate for which the distance to the triangle must be determined.
   * \param[out] dist2Tria Distance squared from the coordinate to the triangle.
   * \param[out] r         Parametric coordinate of the projection.
   * \param[out] s         Parametric coordinate of the projection.
   * \return     True if the projection is inside the triangle and false otherwise.
   */
  bool Dist2ToTriangle(const unsigned long i0, const unsigned long i1, const unsigned long i2, const su2double* coor,
                       su2double& dist2Tria, su2double& r, su2double& s) const;
  /*!
   * \brief Default constructor of the class, disabled.
   */
  CADTElemClass() = delete;
};
