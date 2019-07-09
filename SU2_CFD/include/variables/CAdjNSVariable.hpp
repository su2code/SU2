/*!
 * \file CAdjNSVariable.hpp
 * \brief Main class for defining the variables of the adjoint Navier-Stokes solver.
 * \author F. Palacios, T. Economon
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

#include "CAdjEulerVariable.hpp"

/*!
 * \class CAdjNSVariable
 * \brief Main class for defining the variables of the adjoint Navier-Stokes solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios
 */
class CAdjNSVariable : public CAdjEulerVariable {
private:

public:

  /*!
   * \brief Constructor of the class.
   */
  CAdjNSVariable(void);

  /*!
   * \overload
   * \param[in] val_psirho - Value of the adjoint density (initialization value).
   * \param[in] val_phi - Value of the adjoint velocity (initialization value).
   * \param[in] val_psie - Value of the adjoint energy (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAdjNSVariable(su2double val_psirho, su2double *val_phi, su2double val_psie, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \overload
   * \param[in] val_solution - Pointer to the adjoint value (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAdjNSVariable(su2double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAdjNSVariable(void);

  /*!
   * \brief Set the value of the adjoint velocity.
   * \param[in] val_phi - Value of the adjoint velocity.
   */
  inline void SetPhi_Old(su2double *val_phi) {for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1] = val_phi[iDim]; };

  /*!
   * \brief Set the value of the force projection vector.
   * \param[in] val_ForceProj_Vector - Pointer to the force projection vector.
   */
  inline void SetForceProj_Vector(su2double *val_ForceProj_Vector) {for (unsigned short iDim = 0; iDim < nDim; iDim++) ForceProj_Vector[iDim] = val_ForceProj_Vector[iDim]; }

  /*!
   * \brief Get the value of the force projection vector.
   * \return Pointer to the force projection vector.
   */
  inline su2double *GetForceProj_Vector(void) {return ForceProj_Vector; }

  /*!
   * \brief Set the value of the force projection vector on the solution vector.
   */
  inline void SetVelSolutionOldDVector(void) {for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1] = ForceProj_Vector[iDim]; };

  /*!
   * \brief Set the value of the force projection vector on the old solution vector.
   */
  inline void SetVelSolutionDVector(void) {for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution[iDim+1] = ForceProj_Vector[iDim]; };

};
