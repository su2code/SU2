/*!
 * \file CMeshBoundVariable.hpp
 * \brief Declaration and inlines of the class
 *        to define the variables of the mesh movement at the moving boundaries.
 * \author Ruben Sanchez
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

#include "CMeshVariable.hpp"

class CMeshBoundVariable : public CMeshVariable {
protected:

  su2double *Boundary_Displacement;  /*!< \brief Store the reference coordinates of the mesh. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_coor - Values of the coordinates (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CMeshBoundVariable(su2double *val_coor, unsigned short val_nDim, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CMeshBoundVariable(void);

  /*!
   * \brief Get the value of the displacement imposed at the boundary.
   * \return Value of the boundary displacement.
   */
  inline su2double GetBound_Disp(unsigned short iDim) const final { return Boundary_Displacement[iDim]; }

  /*!
   * \brief Set the boundary displacements.
   * \param[in] val_BoundDisp - Pointer to the boundary displacements.
   */
  inline void SetBound_Disp(const su2double *val_BoundDisp) final {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) Boundary_Displacement[iDim] = val_BoundDisp[iDim];
  }

  /*!
   * \brief Set the boundary displacement.
   * \param[in] iDim - Index of the dimension of interest.
   * \param[in] val_BoundDisp - Value of the boundary displacements.
   */
  inline void SetBound_Disp(unsigned short iDim, const su2double val_BoundDisp) final {
    Boundary_Displacement[iDim] = val_BoundDisp;
  }

  /*!
   * \brief Determine whether the node is a moving vertex.
   * \return True. The node is at the boundary.
   */
  inline bool Get_isVertex(void) const final { return true; }

  /*!
   * \brief Register the boundary displacements of the mesh.
   * \param[in] input - Defines whether we are registering the variable as input or as output.
   */
  inline void Register_BoundDisp(bool input) final {
    if (input) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++)
        AD::RegisterInput(Boundary_Displacement[iVar]);
    }
    else { for (unsigned short iVar = 0; iVar < nVar; iVar++)
        AD::RegisterOutput(Boundary_Displacement[iVar]);
    }
  }

  /*!
   * \brief Recover the value of the adjoint of the boundary displacements.
   */
  inline void GetAdjoint_BoundDisp(su2double *adj_disp) final {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        adj_disp[iVar] = SU2_TYPE::GetDerivative(Boundary_Displacement[iVar]);
    }
  }

};
