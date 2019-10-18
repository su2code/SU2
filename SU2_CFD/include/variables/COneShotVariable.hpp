/*!
 * \file COneShotVariable.hpp
 * \brief Main class for defining the variables of the one-shot optimization solver.
 * \author B. Munguía
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

#include "CDiscAdjVariable.hpp"

/*!
 * \class COneShotVariable
 * \brief Main class for defining the variables of the one-shot optimization solver.
 * \ingroup One_Shot
 * \author B. Munguía.
 */
class COneShotVariable final : public CDiscAdjVariable {
private:
  MatrixType Solution_Store;        /*!< \brief Pointer to store solution of the problem (step k). */
  MatrixType Solution_Former;       /*!< \brief Former solution of the problem (step k-1). */
  MatrixType Solution_Delta;        /*!< \brief Difference of new solution and solution (step k+1 - step k). */
  MatrixType Solution_Delta_Store;  /*!< \brief Difference of solution and old solution (step k - step k-1). */
  MatrixType Solution_Save;         /*!< \brief Pointer to store new solution of the problem (step k+1). */

  MatrixType Sensitivity_ShiftedLagrangian;    /*!< \brief Vector holding the sensitivity of the shifted Lagrangian to the coordinates at this node. */
  MatrixType Sensitivity_AugmentedLagrangian;  /*!< \brief Vector holding the sensitivity of the augmented Lagrangian to the coordinates at this node. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] sol - Pointer to the adjoint value (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  COneShotVariable(const su2double* sol, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~COneShotVariable() = default;

  /*!
   * \brief Set to zero the solution.
   */
  inline void SetSolutionZero(unsigned long iPoint) override {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution(iPoint,iVar) = 0.0;
  }

  /*!
   * \brief Set to zero a particular solution.
   */
  inline void SetSolutionZero(unsigned long iPoint, unsigned short val_var) override {Solution(iPoint,val_var) = 0.0;}

  /*!
   * \brief Set stored variables to the solution.
   */
  inline void Set_StoreSolution(unsigned long iPoint) override {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Store(iPoint,iVar) = Solution(iPoint,iVar);
  }

  /*!
   * \brief Set saved variables to the solution.
   */
  inline void Set_SaveSolution(unsigned long iPoint) override {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Save(iPoint,iVar) = Solution(iPoint,iVar);
  }

  /*!
   * \brief Set former variables to the solution.
   */
  inline void Set_FormerSolution(unsigned long iPoint) override {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Former(iPoint,iVar) = Solution(iPoint,iVar);
  }

  /*!
   * \brief Set former variables to the old solution.
   */
  inline void SetSolution_Former(unsigned long iPoint, su2double *val_solution_old) override {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Former(iPoint,iVar) = val_solution_old[iVar];
  }

  /*!
   * \brief Set stored variables to the old solution.
   */
  inline void SetSolution_Store(unsigned long iPoint, su2double *val_solution_old) override {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Store(iPoint,iVar) = val_solution_old[iVar];
  }

  inline su2double* GetSolution_Store(unsigned long iPoint) override { return Solution_Store[iPoint]; }

  inline su2double GetSolution_Store(unsigned long iPoint, unsigned short iVar) override { return Solution_Store(iPoint,iVar); }

  inline su2double* GetSolution_Save(unsigned long iPoint) override { return Solution_Save[iPoint]; }

  inline su2double GetSolution_Save(unsigned long iPoint, unsigned short iVar) override { return Solution_Save(iPoint,iVar); }

  inline su2double* GetSolution_Delta(unsigned long iPoint) override { return Solution_Delta[iPoint]; }

  inline void SetSolution_Delta(unsigned long iPoint, unsigned short iVar, su2double val_solution_delta) override { Solution_Delta(iPoint,iVar) = val_solution_delta; }

  inline su2double GetSolution_Delta(unsigned long iPoint, unsigned short iVar) override { return Solution_Delta(iPoint,iVar); }

  inline su2double* GetSolution_Former(unsigned long iPoint) override { return Solution_Former[iPoint]; }

  inline su2double GetSolution_Former(unsigned long iPoint, unsigned short iVar) override { return Solution_Former(iPoint,iVar); }

  inline void SetSolution_Delta_Store(unsigned long iPoint, unsigned short iVar, su2double val_solution_delta) override { Solution_Delta_Store(iPoint,iVar) = val_solution_delta; }

  inline su2double GetSolution_Delta_Store(unsigned long iPoint, unsigned short iVar) override { return Solution_Delta_Store(iPoint,iVar); }

  /*!
   * \brief Set the sensitivity in the shifted Lagrangian at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  inline void SetSensitivity_ShiftedLagrangian(unsigned long iPoint, unsigned short iDim, su2double val) override {Sensitivity_ShiftedLagrangian(iPoint,iDim) = val;}

  /*!
   * \brief Set the sensitivity in the augmented Lagrangian at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  inline void SetSensitivity_AugmentedLagrangian(unsigned long iPoint, unsigned short iDim, su2double val) override {Sensitivity_AugmentedLagrangian(iPoint,iDim) = val;}

  /*!
   * \brief Get the Sensitivity in the shifted Lagrangian at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  inline su2double GetSensitivity_ShiftedLagrangian(unsigned long iPoint, unsigned short iDim) override { return Sensitivity_ShiftedLagrangian(iPoint,iDim);}

  /*!
   * \brief Get the Sensitivity in the augmented Lagrangian at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  inline su2double GetSensitivity_AugmentedLagrangian(unsigned long iPoint, unsigned short iDim) override { return Sensitivity_AugmentedLagrangian(iPoint,iDim);}

};