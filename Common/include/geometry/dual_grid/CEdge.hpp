/*!
 * \file CEdge.hpp
 * \brief Headers of the main subroutines for doing the complete dual grid structure.
 *        The subroutines and functions are in the <i>CEdge.cpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 7.0.1 "Blackbird"
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

#include "CDualGrid.hpp"

/*!
 * \class CEdge
 * \brief Class for defining an edge.
 * \author F. Palacios
 */
class CEdge final : public CDualGrid {
private:
  su2double *Coord_CG;      /*!< \brief Center-of-gravity of the element. */
  unsigned long *Nodes;   /*!< \brief Vector to store the global nodes of an element. */
  su2double *Normal;        /*!< \brief Normal al elemento y coordenadas de su centro de gravedad. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_iPoint - First node of the edge.
   * \param[in] val_jPoint - Second node of the edge.
   * \param[in] val_nDim - Number of dimensions of the problem.
   */
  CEdge(unsigned long val_iPoint, unsigned long val_jPoint, unsigned short val_nDim);

  /*!
   * \brief Destructor of the class.
   */
  ~CEdge(void) override;

  /*!
   * \brief Set the center of gravity of the edge.
   * \param[in] val_coord - Coordinates of all the nodes needed for computing the centre of gravity of an edge.
   */
  void SetCoord_CG(su2double **val_coord);

  /*!
   * \brief Obtain the centre of gravity of the edge.
   * \param[in] val_dim - Position to read the coordinate.
   * \return Coordinate <i>val_dim</i> of the centre of gravity.
   */
  inline su2double GetCG(unsigned short val_dim) const { return Coord_CG[val_dim]; }

  /*!
   * \brief Get the nodes of the edge.
   * \param[in] val_node - Position of the node that makes the edge.
   * \return Index of the node that compose the edge.
   */

  inline unsigned long GetNode(unsigned short val_node) const { return Nodes[val_node]; }

  /*!
   * \brief Get the number of nodes of an element.
   * \return Number of nodes that set an edge (2).
   */
  inline unsigned short GetnNodes() const override { return 2; }

  /*!
   * \brief Compute Volume associated to each edge.
   * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
   * \param[in] val_coord_FaceElem_CG - Coordinates of the centre of gravity of the face of an element.
   * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] val_coord_Point - Coordinates of the point that form the control volume.
   * \return Local volume associated to the edge.
   */
  su2double GetVolume(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG, su2double *val_coord_Point) const;

  /*!
   * \overload
   * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
   * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] val_coord_Point - Coordinates of the point that form the control volume.
   * \return Local volume associated to the edge.
   */
  su2double GetVolume(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG, su2double *val_coord_Point) const;

  /*!
   * \brief Set the face that correspond to an edge.
   * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
   * \param[in] val_coord_FaceElem_CG - Coordinates of the centre of gravity of the face of an element.
   * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] config - Definition of the particular problem.
   * \return Compute the normal (dimensional) to the face that makes the control volume boundaries.
   */
  void SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG) override;

  /*!
   * \overload
   * \brief Set the face that correspond to an edge.
   * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
   * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] config - Definition of the particular problem.
   * \return Compute the normal (dimensional) to the face that makes the contorl volume boundaries.
   */
  void SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG) override;

  /*!
   * \brief Copy the the normal vector of a face.
   * \param[in] val_normal - Vector where the subroutine is goint to copy the normal (dimensional).
   */
  inline void GetNormal(su2double *val_normal) const override {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      val_normal[iDim] = Normal[iDim];
  }

  /*!
   * \brief Get the normal to a face of the control volume asociated with an edge.
   * \return Dimensional normal vector, the modulus is the area of the face.
   */
  inline su2double *GetNormal(void) override {  return Normal; }

  /*!
   * \brief Initialize normal vector.
   */
  inline void SetZeroValues(void) override {
    for (unsigned short iDim = 0; iDim < nDim; iDim ++)
      Normal[iDim] = 0.0;
  }

  /*!
   * \brief Set the normal vector.
   * \param[in] val_face_normal - Vector to initialize the normal vector.
   * \return Value of the normal vector.
   */
  inline void SetNormal(const su2double *val_face_normal) override {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Normal[iDim]=val_face_normal[iDim];
  }

  /*!
   * \brief Add a vector to the normal vector.
   * \param[in] val_face_normal - Vector to add to the normal vector.
   */
  inline void AddNormal(const su2double *val_face_normal) override {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Normal[iDim] += val_face_normal[iDim];
  }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline su2double *GetCoord(void) override { return NULL; }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline void SetCoord(const su2double *val_coord) override { }

};
