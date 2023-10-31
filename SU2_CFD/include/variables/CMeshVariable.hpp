/*!
 * \file CMeshVariable.hpp
 * \brief Declaration and inlines of the class
 *        to define the variables of the mesh movement.
 * \author Ruben Sanchez
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

#include "CVariable.hpp"

class CMeshVariable : public CVariable {
protected:

  VectorType WallDistance;  /*!< \brief Store the wall distance in reference coordinates. */
  MatrixType Mesh_Coord;    /*!< \brief Store the reference coordinates of the mesh. */

  /*!
   * \brief Constructor of the class.
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CMeshVariable(unsigned long npoint, unsigned long ndim, CConfig *config);

public:
  /*!
   * \brief Destructor of the class.
   */
  ~CMeshVariable() override = default;

  /*!
   * \brief Get the value of the undeformed coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \return Value of the original coordinate iDim.
   */
  inline su2double GetMesh_Coord(unsigned long iPoint, unsigned long iDim) const final { return Mesh_Coord(iPoint,iDim); }

  /*!
   * \brief Get the undeformed coordinates.
   * \return Pointer to the reference coordinates.
   */
  inline const su2double *GetMesh_Coord(unsigned long iPoint) const final { return Mesh_Coord[iPoint]; }

  /*!
   * \brief Get the undeformed coordinates for the entire domain.
   */
  inline const MatrixType* GetMesh_Coord() const final { return &Mesh_Coord; }

  /*!
   * \brief Set the value of the undeformed coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \param[in] val_coord - Value of Mesh_Coord[nDim]
   */
  inline void SetMesh_Coord(unsigned long iPoint, unsigned long iDim, su2double val_coord) final {
    Mesh_Coord(iPoint,iDim) = val_coord;
  }

  /*!
   * \brief Get the value of the wall distance in reference coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \return Value of the wall distance in reference coordinates.
   */
  inline su2double GetWallDistance(unsigned long iPoint) const final { return WallDistance(iPoint); }

  /*!
   * \brief Set the value of the wall distance in reference coordinates.
   * \param[in] val_dist - Value of wall distance.
   */
  inline void SetWallDistance(unsigned long iPoint, su2double val_dist) final { WallDistance(iPoint) = val_dist; }

  /*!
   * \brief Register the reference coordinates of the mesh.
   */
  void Register_MeshCoord() final;

  /*!
   * \brief Recover the value of the adjoint of the mesh coordinates.
   */
  inline void GetAdjoint_MeshCoord(unsigned long iPoint, su2double *adj_mesh) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) {
      adj_mesh[iDim] = SU2_TYPE::GetDerivative(Mesh_Coord(iPoint,iDim));
      AD::ResetInput(Mesh_Coord(iPoint,iDim));
    }
  }

};
