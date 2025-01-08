/*!
 * \file CTurboVertex.hpp
 * \brief Headers of the main subroutines for doing the complete dual grid structure.
 *        The subroutines and functions are in the <i>CTurboVertex.cpp</i> file.
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

#include "CVertex.hpp"

/*!
 * \class CTurboVertex
 * \brief Class for vertex definition for turbomachinery (equivalent to edges, but for the boundaries).
 * \author S. Vitale
 */
class CTurboVertex final : public CVertex {
 private:
  su2double* TurboNormal; /*!< \brief Normal for computing correct turbomachinery quantities. */
  su2double Area;         /*!< \brief Value of the face area associated to the vertex */
  //  su2double PitchCoord;       /*!< \brief Value of the abscissa pitch wise */
  su2double AngularCoord;      /*!< \brief Value of the angular coordinate  */
  su2double DeltaAngularCoord; /*!< \brief Value of the angular coordinate w.r.t. the minimum pitch point  */
  su2double RelAngularCoord;   /*!< \brief Value of the angular coordinate w.r.t. the minimum pitch point  */

  unsigned long OldVertex; /*!< \brief Value of the vertex numeration before the ordering */
  int GlobalIndex;         /*!< \brief Value of the vertex numeration after the ordering and global with respect to MPI
                              partinioning */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_point - Node of the vertex.
   * \param[in] val_nDim - Number of dimensions of the problem.
   */
  CTurboVertex(unsigned long val_point, unsigned short val_nDim);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurboVertex(void) override;

  /*!
   * \brief set Normal in the turbomachinery frame of reference.
   * \param[in] val_normal - normal vector.
   */
  inline void SetTurboNormal(const su2double* val_normal) {
    unsigned short iDim;
    for (iDim = 0; iDim < nDim; iDim++) TurboNormal[iDim] = val_normal[iDim];
  }

  /*!
   * \brief set face Area.
   * \param[in] val_area - value of the face area.
   */
  inline void SetArea(su2double val_area) { Area = val_area; }

  /*!
   * \brief get face Area associate to the vertex.
   */
  inline su2double GetArea(void) const { return Area; }

  /*!
   * \brief Copy the the turbo normal vector of a face.
   * \param[in] val_normal - Vector where the subroutine is goint to copy the normal (dimensionaless).
   */
  inline void GetTurboNormal(su2double* val_normal) const {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) val_normal[iDim] = TurboNormal[iDim];
  }

  /*!
   * \brief Get the turbo normal to a face where turboperformance are computed .
   * \return Dimensionaless normal vector, the modulus is the area of the face.
   */
  inline su2double* GetTurboNormal(void) { return TurboNormal; }

  /*!
   * \brief set vertex value not ordered.
   * \param[in] val_vertex - value of the vertex before ordering.
   */
  inline void SetOldVertex(unsigned long val_vertex) { OldVertex = val_vertex; }

  /*!
   * \brief retrieve vertex value not ordered.
   */
  inline unsigned long GetOldVertex(void) const { return OldVertex; }

  /*!
   * \brief set global index for ordered span-wise turbovertex.
   */
  inline void SetGlobalVertexIndex(int globalindex) { GlobalIndex = globalindex; }

  /*!
   * \brief get global index for ordered span-wise turbovertex.
   */
  inline int GetGlobalVertexIndex(void) const { return GlobalIndex; }

  /*!
   * \brief set angular coord.
   */
  inline void SetAngularCoord(su2double angCoord) { AngularCoord = angCoord; }

  /*!
   * \brief get angular coord.
   */
  inline su2double GetAngularCoord(void) const { return AngularCoord; }

  /*!
   * \brief set angular coord.
   */
  inline void SetDeltaAngularCoord(su2double deltaAngCoord) { DeltaAngularCoord = deltaAngCoord; }

  /*!
   * \brief get angular coord.
   */
  inline su2double GetDeltaAngularCoord(void) const { return DeltaAngularCoord; }

  /*!
   * \brief set angular coord.
   */
  inline void SetRelAngularCoord(su2double minAngCoord) { RelAngularCoord = AngularCoord - minAngCoord; }

  /*!
   * \brief get angular coord.
   */
  inline su2double GetRelAngularCoord(void) const { return RelAngularCoord; }
};
