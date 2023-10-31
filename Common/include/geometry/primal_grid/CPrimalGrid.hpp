/*!
 * \file CPrimalGrid.hpp
 * \brief Headers of the main subroutines for storing the primal grid structure.
 *        The subroutines and functions are in the <i>primal_grid_structure.cpp</i> file.
 * \author F. Palacios
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

#include <iostream>
#include <vector>
#include <limits>
#include <cstdlib>
#include <limits>
#include <memory>

#include "../../option_structure.hpp"

/*!
 * \class CPrimalGrid
 * \brief Class to define the numerical primal grid.
 * \author F. Palacios, T. Economon, M. Aehle.
 */
class CPrimalGrid {
 protected:
  /* If this is a domain element, stores the global index.
   * If this is a boundary element, stores the index of the adjacent domain element. */
  unsigned long GlobalIndex_DomainElement;

  std::unique_ptr<unsigned long[]> Nodes;    /*!< \brief Global node indices of the element. */
  std::unique_ptr<long[]> Neighbor_Elements; /*!< \brief Vector to store the elements surronding this element. */

  su2double Coord_CG[3] = {0.0}; /*!< \brief Coordinates of the center-of-gravity of the element. */
  su2double Volume;              /*!< \brief Volume of the element. */
  su2double LenScale;            /*!< \brief Length scale of the element. */

  unsigned short TimeLevel; /*!< \brief Time level of the element for time accurate local time stepping. */
  /*!< \brief Vector to store the periodic index of a neighbor, -1 indicates no periodic transformation to the neighbor.
   */
  int8_t PeriodIndexNeighbors[N_FACES_MAXIMUM];

  /*! \brief Whether or not the Jacobian of the faces can be considered
   * constant in the transformation to the standard element. */
  bool JacobianFaceIsConstant[N_FACES_MAXIMUM];
  bool ElementOwnsFace[N_FACES_MAXIMUM]; /*!< \brief Whether or not the element owns each face. */
  const bool FEM;                        /*!< \brief Whether this is a FEM element. */

 public:
  CPrimalGrid() = delete;

  /*!
   * \brief Constructor of the class.
   * \param[in] FEM - Whether this is a FEM element.
   * \param[in] nNodes - Number of nodes.
   * \param[in] nNeighbor_Elements - Number of neighbor elements.
   */
  CPrimalGrid(bool FEM, unsigned short nNodes, unsigned short nNeighbor_Elements);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CPrimalGrid() = default;

  /*!
   * \brief Get the nodes shared by the primal grid element.
   * \param[in] val_node - Local (to the element) index of the node (lines have 2 nodes, triangles 3 nodes etc).
   * \return Global index of the node.
   */
  inline unsigned long GetNode(unsigned short val_node) const { return Nodes[val_node]; }

  /*!
   * \brief Set the nodes shared by the primal grid element.
   * \param[in] val_node - Local (to the element) index of the node (lines have 2 nodes, triangles 3 nodes etc).
   * \param[in] val_point - Global index of the node.
   */
  inline void SetNode(unsigned short val_node, unsigned long val_point) {
    assert(val_node < GetnNodes());
    Nodes[val_node] = val_point;
  }

  /*!
   * \brief Get the elements that surround an element.
   * \param[in] val_face - Local index of the face.
   * \return Global index of the element.
   */
  inline long GetNeighbor_Elements(unsigned short val_face) const { return Neighbor_Elements[val_face]; }

  /*!
   * \brief Set the elements that surround an element.
   * \param[in] val_elem - Global index of the element.
   * \param[in] val_face - Local index of the face.
   */
  inline void SetNeighbor_Elements(unsigned long val_elem, unsigned short val_face) {
    Neighbor_Elements[val_face] = val_elem;
  }

  /*!
   * \brief Make available the length scale of the element.
   * \return The length scale of the element.
   */
  inline su2double GetLengthScale() const { return LenScale; }

  /*!
   * \brief Set the length scale of the element.
   * \param[in] val_lenScale - Length scale of the element.
   */
  inline void SetLengthScale(su2double val_lenScale) { LenScale = val_lenScale; }

  /*!
   * \brief Make available the time level of the element.
   * \return The time level of the element.
   */
  inline unsigned short GetTimeLevel() const { return TimeLevel; }

  /*!
   * \brief Set the time level of the element.
   * \param[in] val_timeLevel - Time level of the element.
   */
  inline void SetTimeLevel(unsigned short val_timeLevel) { TimeLevel = val_timeLevel; }

  /*!
   * \brief Get the boolean to indicate whether or not this element owns the face
   *        between the current and the adjacent element with index val_face.
   * \param[in] val_face - Local index of the face.
   * \return   Boolean to indicate whether or not the face is owned by this element.
   */
  inline bool GetOwnerFace(unsigned short val_face) const { return ElementOwnsFace[val_face]; }

  /*!
   * \brief Set the boolean to indicate whether or not this element owns the face
   *        between the current and the adjacent element with index val_face.
   * \param[in] val_owner - Whether or not this element owns the face.
   * \param[in] val_face  - Local index of the face.
   */
  inline void SetOwnerFace(bool val_owner, unsigned short val_face) { ElementOwnsFace[val_face] = val_owner; }

  /*!
   * \brief Get the index of the periodic transformation to the neighboring element.
   * \param[in] val_face - Local index of the face.
   * \return   Index of the periodic transformation to the neighboring element.
   */
  inline short GetPeriodicIndex(unsigned short val_face) const { return PeriodIndexNeighbors[val_face]; }

  /*!
   * \brief Set the index of the periodic transformation to the neighboring element.
   * \param[in] val_periodic - Index of the periodic marker to which the face belongs.
   * \param[in] val_face     - Local index of the face.
   */
  inline void SetPeriodicIndex(unsigned short val_periodic, unsigned short val_face) {
    assert(val_periodic < std::numeric_limits<int8_t>::max());
    PeriodIndexNeighbors[val_face] = val_periodic;
  }

  /*!
   * \brief Get whether or not the Jacobian of the given face is considered constant.
   * \param[in] val_face - Local index of the face.
   * \return  Whether or not the Jacobian of the face is considered constant.
   */
  inline bool GetJacobianConstantFace(unsigned short val_face) const { return JacobianFaceIsConstant[val_face]; }

  /*!
   * \brief Set whether or not the Jacobian of the given face is considered constant.
   * \param[in] val_JacFaceIsConstant - Boolean to indicate whether or not the Jacobian is constant.
   * \param[in] val_face              - Local index of the face.
   */
  inline void SetJacobianConstantFace(bool val_JacFaceIsConstant, unsigned short val_face) {
    JacobianFaceIsConstant[val_face] = val_JacFaceIsConstant;
  }

  /*!
   * \brief Set the center of gravity of an element (including edges).
   * \param[in] nDim - Number of dimensions (2 or 3).
   * \param[in] val_coord - Coordinates of the element.
   */
  template <class T>
  inline su2double* SetCoord_CG(unsigned short nDim, const T& val_coord) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Coord_CG[iDim] = 0.0;
      for (unsigned short iNode = 0; iNode < GetnNodes(); iNode++)
        Coord_CG[iDim] += val_coord[iNode][iDim] / su2double(GetnNodes());
    }
    return Coord_CG;
  }

  /*!
   * \brief Get the center of gravity of an element (including edges).
   * \param[in] val_dim - Coordinate of the center of gravity.
   * \return Coordinates of the center of gravity.
   */
  inline su2double GetCG(unsigned short val_dim) const { return Coord_CG[val_dim]; }
  inline const su2double* GetCG() const { return Coord_CG; }

  /*!
   * \brief Set the center of gravity of an element (including edges).
   * \param[in] val_coord - Coordinates of the element.
   */
  inline void SetVolume(su2double val_volume) { Volume = val_volume; }

  /*!
   * \brief Get the center of gravity of an element (including edges).
   * \param[in] val_dim - Coordinate of the center of gravity.
   * \return Coordinates of the center of gravity.
   */
  inline su2double GetVolume() const { return Volume; }

  /*!
   * \brief Reset the array, which stores whether or not the faces have a constant Jacobian, to false.
   */
  inline void ResetJacobianConstantFaces() {
    for (auto i = 0; i < N_FACES_MAXIMUM; ++i) JacobianFaceIsConstant[i] = false;
  }

  /*!
   * \brief Initialize the information about the neighboring elements.
   * \param[in] val_nFaces - Number of faces for which neighboring information must be initialized.
   */
  void InitializeNeighbors(unsigned short val_nFaces);

  /*!
   * \brief A virtual member.
   * \param[in] val_color - New color of the element.
   */
  inline virtual void SetColor(unsigned long val_color) {}

  /*!
   * \brief A virtual member.
   * \return The color of the element in the partitioning.
   */
  inline virtual unsigned long GetColor() const { return std::numeric_limits<unsigned long>::max(); }

  /*!
   * \brief Get the element global index in a parallel computation.
   * \return Global index of the element in a parallel computation.
   */
  inline unsigned long GetGlobalIndex() const { return GlobalIndex_DomainElement; }

  /*!
   * \brief Set the global index for an element in a parallel computation.
   * \return Global index of an element in a parallel computation.
   */
  inline void SetGlobalIndex(unsigned long val_globalindex) { GlobalIndex_DomainElement = val_globalindex; }

  /*!
   * \brief Set the index of the domain element of which this boundary element is a face.
   * \param[in] val_domainelement - Value to set.
   */
  inline void SetDomainElement(unsigned long val_domainelement) { GlobalIndex_DomainElement = val_domainelement; }

  /*!
   * \brief Get the index of the domain element of which this boundary element is a face.
   */
  inline unsigned long GetDomainElement() const { return GlobalIndex_DomainElement; }

  /*!
   * \brief A pure virtual member.
   */
  virtual void Change_Orientation() = 0;

  /*!
   * \brief A pure virtual member.
   * \return Type of the element using VTK nomenclature.
   */
  inline virtual unsigned short GetRotation_Type() const { return 0; }

  /*!
   * \brief A pure virtual member.
   * \param[in] val_rotation_type - Kind of rotation/traslation that must be applied.
   */
  inline virtual void SetRotation_Type(unsigned short val_rotation_type) {}

  /*-- The following pure virtual functions are overridden in
   * CPrimalGridWithConnectivity, except for the FEM classes. --*/

  /*!
   * \brief Get number of nodes of the element.
   * \return Number of nodes.
   */
  virtual unsigned short GetnNodes() const = 0;

  /*!
   * \brief Get number of faces of the element.
   * \return Number of faces.
   */
  virtual unsigned short GetnFaces() const = 0;

  /*!
   * \brief Get number of nodes of a face of the element.
   * \param[in] val_face - Index of the face among all faces of the element.
   * \return Number of nodes contained in the face.
   */
  inline virtual unsigned short GetnNodesFace(unsigned short val_face) const { return 0; }

  /*!
   * \brief Get the maximum number of nodes contained in a face of the element.
   * \return Maximum number of nodes contained in a face.
   */
  virtual unsigned short GetMaxNodesFace() const = 0;

  /*!
   * \brief Get nodes contained in a face.
   * \param[in] val_face - Index of the face among all faces of the element.
   * \param[in] val_index - Index of the node among all nodes of the face.
   * \return - Index of the node among all nodes of the element.
   */
  virtual unsigned short GetFaces(unsigned short val_face, unsigned short val_index) const = 0;

  /*!
   * \brief Get number of neighbor nodes of a node.
   * \param[in] val_node - Index of node among all nodes of the element.
   * \return Number of neighbor nodes.
   */
  virtual unsigned short GetnNeighbor_Nodes(unsigned short val_node) const = 0;

  /*!
   * \brief Get neighbor nodes of a node.
   * \param[in] val_node - Index of node N among all nodes of the element.
   * \param[in] val_index - Index of the neighbor node among all neighbor nodes of the node N.
   * \return Index of neighbor node among all nodes of the element.
   */
  virtual unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) const = 0;

  /*!
   * \brief Get VTK type.
   * \return Type of element using the VTK nomenclature.
   */
  virtual unsigned short GetVTK_Type() const = 0;

  /*-- Until here --*/

  /*!
   * \brief Virtual function, that must be overwritten by the derived class, if needed.
   * \param[out] nFaces         - Number of faces of this element.
   * \param[out] nPointsPerFace - Number of corner points for each of the faces.
   * \param[out] faceConn       - Global IDs of the corner points of the faces.
   */
  inline virtual void GetCornerPointsAllFaces(unsigned short& nFaces, unsigned short nPointsPerFace[],
                                              unsigned long faceConn[6][4]) const {}

  /*!
   * \brief Virtual function to make available the global ID of this element.
   * \return The global ID of this element.
   */
  inline virtual unsigned long GetGlobalElemID() const { return 0; }

  /*!
   * \brief Virtual function to make available the global offset of the solution DOFs.
   * \return The global offset of the solution DOFs.
   */
  inline virtual unsigned long GetGlobalOffsetDOFsSol() const { return 0; }

  /*!
   * \brief Virtual function to make available the polynomial degree of the grid.
   * \return The polynomial degree of the grid.
   */
  inline virtual unsigned short GetNPolyGrid() const { return 0; }

  /*!
   * \brief Virtual function to make available the polynomial degree of the solution.
   * \return The polynomial degree of the solution.
   */
  inline virtual unsigned short GetNPolySol() const { return 0; }

  /*!
   * \brief Virtual function to make available the number of DOFs of the grid in the element.
   * \return The number of DOFs of the Grid in the element.
   */
  inline virtual unsigned short GetNDOFsGrid() const { return 0; }

  /*!
   * \brief Virtual function to make available the number of DOFs of the solution in the element.
   * \return The number of DOFs of the solution in the element.
   */
  inline virtual unsigned short GetNDOFsSol() const { return 0; }

  /*!
   * \brief Virtual function to get whether or not the Jacobian is considered constant.
   * \return True if the Jacobian is (almost) constant and false otherwise.
   */
  inline virtual bool GetJacobianConsideredConstant() const { return false; }

  /*!
   * \brief Virtual function to set the value of JacobianConsideredConstant.
   * \param[in] val_JacobianConsideredConstant - The value to be set for JacobianConsideredConstant.
   */
  inline virtual void SetJacobianConsideredConstant(bool val_JacobianConsideredConstant) {}

  /*!
   * \brief Virtual function to correct the offset of the global DOFs.
   * \param[in] val_offsetRank - The offset that must be added for this rank.
   */
  inline virtual void AddOffsetGlobalDOFs(const unsigned long val_offsetRank) {}

  /*!
   * \brief Virtual function to add the given donor ID to the donor elements for the wall function treatment.
   * \param[in] donorElement - Element to be added to donor elements.
   */
  inline virtual void AddDonorWallFunctions(const unsigned long donorElement) {}

  /*!
   * \brief Virtual function to make available the number of donor elements for the wall function treatment.
   * \return The number of donor elements.
   */
  inline virtual unsigned short GetNDonorsWallFunctions() const { return 0; }

  /*!
   * \brief Virtual function to make available the pointer to the vector for the donor elements
            for the wall function treatment.
   * \return The pointer to the data of donorElementsWallFunctions.
   */
  inline virtual unsigned long* GetDonorsWallFunctions() { return nullptr; }
  inline virtual const unsigned long* GetDonorsWallFunctions() const { return nullptr; }

  /*!
   * \brief Virtual function to set the global ID's of the donor elements for the wall function treatment.
   * \param[in] donorElements - Vector, which contain the donor elements.
   */
  inline virtual void SetDonorsWallFunctions(const std::vector<unsigned long>& donorElements) {}

  /*!
   * \brief Virtual function to remove the multiple donors for the wall function treatment.
   */
  inline virtual void RemoveMultipleDonorsWallFunctions() {}
};

/*! \class CPrimalGridWithConnectivity
 * \brief Override the connectivity getters of CPrimalGrid.
 *
 * \details Non-FEM primal grid classes like CLine, CTriangle etc
 * are derived from CPrimalGridWithConnectivity<CLineConnectivity> etc.
 * The Connectivity class must have the following static members of
 * type unsigned short / unsigned short array, typically constexpr:
 * - nNodes: Number of nodes of the element.
 * - nFaces:  Number of faces of the element.
 * - nNodesFace[]: Number of nodes of each face of the element.
 * - maxNodesFace: Maximum number of nodes for a face.
 * - Faces[][]: Matrix to store the local nodes of all the faces.
 * - nNeighbor_Nodes[]: Number of neighbor nodes of each node of the element.
 * - Neighbor_Nodes[][]: Matrix to store the neighbors of all the nodes.
 * - VTK_Type: Type of element using VTK nomenclature.
 *
 * The getter functions for connectivity in CPrimalGrid have final overrides
 * in this class, accessing the static members of Connectivity.
 *
 * \tparam Connectivity - class defining the connectivity structure
 */
template <typename Connectivity>
class CPrimalGridWithConnectivity : public CPrimalGrid {
 public:
  CPrimalGridWithConnectivity(bool FEM) : CPrimalGrid(FEM, Connectivity::nNodes, Connectivity::nFaces) {}

  inline unsigned short GetnNodes() const final { return Connectivity::nNodes; }

  inline unsigned short GetnFaces() const final { return Connectivity::nFaces; }

  inline unsigned short GetnNodesFace(unsigned short val_face) const final {
    assert(val_face < Connectivity::nFaces);
    return Connectivity::nNodesFace[val_face];
  }

  inline unsigned short GetMaxNodesFace() const final { return Connectivity::maxNodesFace; }

  inline unsigned short GetFaces(unsigned short val_face, unsigned short val_index) const final {
    assert(val_face < GetnFaces() && val_index < GetnNodesFace(val_face));
    return Connectivity::Faces[val_face][val_index];
  }

  inline unsigned short GetnNeighbor_Nodes(unsigned short val_node) const final {
    assert(val_node < Connectivity::nNodes);
    return Connectivity::nNeighbor_Nodes[val_node];
  }

  inline unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) const final {
    assert(val_node < Connectivity::nNodes && val_index < GetnNeighbor_Nodes(val_node));
    return Connectivity::Neighbor_Nodes[val_node][val_index];
  }

  inline unsigned short GetVTK_Type() const final { return Connectivity::VTK_Type; }
};
