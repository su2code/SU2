/*!
 * \file CVertex.hpp
 * \brief Headers of the main subroutines for doing the complete dual grid structure.
 *        The subroutines and functions are in the <i>CVertex.cpp</i> file.
 * \author F. Palacios, T. Economon
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

#include "CDualGrid.hpp"

/*!
 * \class CVertex
 * \brief Class for vertex definition (equivalent to edges, but for the boundaries).
 * \author F. Palacios
 */
class CVertex : public CDualGrid {
 protected:
  unsigned long Nodes[1];         /*!< \brief Vector to store the global nodes of an element. */
  su2double Normal[3] = {0.0};    /*!< \brief Normal coordinates of the element and its center of gravity. */
  su2double Aux_Var;              /*!< \brief Auxiliar variable defined only on the surface. */
  su2double CartCoord[3] = {0.0}; /*!< \brief Vertex cartesians coordinates. */
  su2double VarCoord[3] = {0.0}; /*!< \brief Used for storing the coordinate variation due to a surface modification. */
  long PeriodicPoint[5] = {-1};  /*!< \brief Store the periodic point of a boundary (iProcessor, iPoint) */
  bool ActDisk_Perimeter = false;      /*!< \brief Identify nodes at the perimeter of the actuator disk */
  short Rotation_Type;                 /*!< \brief Type of rotation associated with the vertex (MPI and periodic) */
  unsigned long Normal_Neighbor;       /*!< \brief Index of the closest neighbor. */
  su2double Basis_Function[3] = {0.0}; /*!< \brief Basis function values for interpolation across zones. */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_point - Node of the vertex.
   * \param[in] val_nDim - Number of dimensions of the problem.
   */
  CVertex(unsigned long val_point, unsigned short val_nDim);

  /*!
   * \brief Get the number of nodes of a vertex.
   * \return Number of nodes that set a vertex (1).
   */
  inline unsigned short GetnNodes() const override { return 1; }

  /*!
   * \brief Get the node of the vertex.
   * \return Index of the node that compose the vertex.
   */
  inline unsigned long GetNode() const { return Nodes[0]; }

  /*!
   * \brief Set the face that correspond to a vertex.
   * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
   * \param[in] val_coord_FaceElem_CG - Coordinates of the centre of gravity of the face of an element.
   * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \return Compute the normal (dimensional) to the face that makes the vertex.
   */
  void SetNodes_Coord(const su2double* val_coord_Edge_CG, const su2double* val_coord_FaceElem_CG,
                      const su2double* val_coord_Elem_CG) override;

  /*!
   * \overload
   * \param[in] val_coord_Edge_CG - Coordinates of the centre of gravity of the edge.
   * \param[in] val_coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \return Compute the normal (dimensional) to the face that makes the vertex.
   */
  void SetNodes_Coord(const su2double* val_coord_Edge_CG, const su2double* val_coord_Elem_CG) override;

  /*!
   * \brief Copy the the normal vector of a face.
   * \param[in] val_normal - Vector where the subroutine is goint to copy the normal (dimensional).
   */
  inline void GetNormal(su2double* val_normal) const override {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) val_normal[iDim] = Normal[iDim];
  }

  /*!
   * \brief Get the normal to a face of the control volume asociated with a vertex.
   * \return Dimensional normal vector, the modulus is the area of the face.
   */
  inline su2double* GetNormal(void) override { return Normal; }

  /*!
   * \brief Get the ith component of the normal.
   */
  inline su2double GetNormal(unsigned short iDim) const { return Normal[iDim]; }

  /*!
   * \brief Initialize normal vector.
   */
  inline void SetZeroValues(void) override {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) Normal[iDim] = 0.0;
  }

  /*!
   * \brief Set the value of an auxiliary variable for gradient computation.
   * \param[in] val_auxvar - Value of the auxiliar variable.
   */
  inline void SetAuxVar(su2double val_auxvar) { Aux_Var = val_auxvar; }

  /*!
   * \brief Get the value of an auxiliary variable for gradient computation.
   * \return Value of the auxiliar variable.
   */
  inline su2double GetAuxVar(void) const { return Aux_Var; }

  /*!
   * \brief Add the value of an auxiliary variable for gradient computation.
   * \param[in] val_auxvar - Value of the auxiliar variable.
   */
  inline void AddAuxVar(su2double val_auxvar) { Aux_Var += val_auxvar; }

  /*!
   * \brief Set the normal vector.
   * \param[in] val_face_normal - Vector to initialize the normal vector.
   * \return Value of the normal vector.
   */
  inline void SetNormal(const su2double* val_face_normal) override {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) Normal[iDim] = val_face_normal[iDim];
  }

  /*!
   * \brief Add a vector to the normal vector.
   * \param[in] val_face_normal - Vector to add to the normal vector.
   */
  inline void AddNormal(const su2double* val_face_normal) override {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) Normal[iDim] += val_face_normal[iDim];
  }

  /*!
   * \brief Set the value of the coordinate variation due to a surface modification.
   * \param[in] val_varcoord - Variation of the coordinate.
   */
  inline void SetVarCoord(const su2double* val_varcoord) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) VarCoord[iDim] = val_varcoord[iDim];
  }

  /*!
   * \brief Add the value of the coordinate variation due to a surface modification.
   * \param[in] val_varcoord - Variation of the coordinate.
   */
  inline void AddVarCoord(const su2double* val_varcoord) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) VarCoord[iDim] += val_varcoord[iDim];
  }

  /*!
   * \brief Get the value of the coordinate variation due to a surface modification.
   * \return Variation of the coordinate.
   */
  inline su2double* GetVarCoord(void) { return VarCoord; }

  /*!
   * \brief Set the value of the cartesian coordinate for the vertex.
   * \param[in] val_coord - Value of the cartesian coordinate.
   */
  inline void SetCoord(const su2double* val_coord) override {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) CartCoord[iDim] = val_coord[iDim];
  }

  /*!
   * \brief Get the value of the cartesian coordinate for the vertex.
   * \return Value of the cartesian coordinate of the vertex.
   */
  inline su2double* GetCoord(void) override { return CartCoord; }

  /*!
   * \brief Get the value of the cartesian coordinate for the vertex.
   * \param[in] val_dim - Variable of the dimension.
   * \return Value of the cartesian coordinate of the vertex.
   */
  inline su2double GetCoord(unsigned short val_dim) const { return CartCoord[val_dim]; }

  /*!
   * \brief Set the type of rotation associated to the vertex.
   * \param[in] val_rotation_type - Value of the rotation that will be applied to the solution at the vertex
   */
  inline void SetRotation_Type(short val_rotation_type) { Rotation_Type = val_rotation_type; }

  /*!
   * \brief Get the type of rotation associated to the vertex.
   * \return Value of the rotation that must be applied to the solution of the vertex
   */
  inline short GetRotation_Type(void) const { return Rotation_Type; }

  /*!
   * \overload
   * \param[in] val_periodicpoint - Value of periodic point of the vertex.
   * \param[in] val_processor - Processor where the point belong.
   */
  inline void SetDonorPoint(long val_periodicpoint, long val_processor) {
    PeriodicPoint[0] = val_periodicpoint;
    PeriodicPoint[1] = val_processor;
    PeriodicPoint[2] = 0;
  }

  /*!
   * \overload
   * \param[in] val_periodicpoint - Value of periodic point of the vertex.
   * \param[in] val_processor - Processor where the point belong.
   */
  inline void SetDonorPoint(long val_periodicpoint, long val_periodicglobalindex, long val_periodicvertex,
                            long val_periodicmarker, long val_processor) {
    PeriodicPoint[0] = val_periodicpoint;
    PeriodicPoint[1] = val_processor;
    PeriodicPoint[2] = val_periodicglobalindex;
    PeriodicPoint[3] = val_periodicvertex;
    PeriodicPoint[4] = val_periodicmarker;
  }

  /*!
   * \overload
   * \param[in] val_periodicpoint - Value of periodic point of the vertex.
   * \param[in] val_processor - Processor where the point belong.
   * \param[in] val_globalindex - Global index of the donor point.
   */
  inline void SetDonorPoint(long val_periodicpoint, long val_processor, long val_globalindex) {
    PeriodicPoint[0] = val_periodicpoint;
    PeriodicPoint[1] = val_processor;
    PeriodicPoint[2] = val_globalindex;
  }

  /*!
   * \overload
   * \param[in] val_periodicpoint - Value of periodic point of the vertex.
   * \param[in] val_processor - Processor where the point belong.
   */
  inline void SetActDisk_Perimeter(bool val_actdisk_perimeter) { ActDisk_Perimeter = val_actdisk_perimeter; }

  /*!
   * \brief Get the value of the periodic point of a vertex.
   * \return Value of the periodic point of a vertex.
   */
  inline long GetDonorPoint(void) const { return PeriodicPoint[0]; }

  /*!
   * \brief Get the value of the periodic point of a vertex.
   * \return Value of the periodic point of a vertex.
   */
  inline long GetDonorMarker(void) const { return PeriodicPoint[4]; }

  /*!
   * \brief Get the value of the periodic point of a vertex.
   * \return Value of the periodic point of a vertex.
   */
  inline long GetDonorVertex(void) const { return PeriodicPoint[3]; }

  /*!
   * \brief Get the value of the periodic point of a vertex.
   * \return Value of the periodic point of a vertex.
   */
  inline long GetDonorGlobalIndex(void) const { return PeriodicPoint[2]; }

  /*!
   * \brief Get the value of the periodic point of a vertex.
   * \return Value of the periodic point of a vertex.
   */
  inline long GetGlobalDonorPoint(void) const { return PeriodicPoint[2]; }

  /*!
   * \brief Get the value of the periodic point of a vertex.
   * \return Value of the periodic point of a vertex.
   */
  inline long GetDonorProcessor(void) const { return PeriodicPoint[1]; }

  /*!
   * \brief Get the value of the periodic point of a vertex, and its somain
   * \return Value of the periodic point of a vertex, and the domain.
   */
  inline long* GetPeriodicPointDomain(void) { return PeriodicPoint; }

  /*!
   * \brief Get the value of the periodic point of a vertex, and its somain
   * \return Value of the periodic point of a vertex, and the domain.
   */
  inline bool GetActDisk_Perimeter(void) const { return ActDisk_Perimeter; }

  /*!
   * \brief Set the finite element basis functions needed for interpolation.
   * \param[in] val_node - a node index of the owner element.
   * \param[in] val_basis - basis function value for the node.
   */
  inline void SetBasisFunction(unsigned short val_node, su2double val_basis) { Basis_Function[val_node] = val_basis; }

  /*!
   * \brief Get the finite element basis functions needed for interpolation.
   * \param[in] val_node - a node index of the owner element.
   * \return Value of the basis function for this node.
   */
  inline su2double GetBasisFunction(unsigned short val_node) { return Basis_Function[val_node]; }

  /*!
   * \brief Set the index of the closest neighbor to a point on the boundaries.
   * \param[in] val_Normal_Neighbor - Index of the closest neighbor.
   */
  inline void SetNormal_Neighbor(unsigned long val_Normal_Neighbor) { Normal_Neighbor = val_Normal_Neighbor; }

  /*!
   * \brief Get the value of the closest neighbor.
   * \return Index of the closest neighbor.
   */
  inline unsigned long GetNormal_Neighbor(void) const { return Normal_Neighbor; }
};
