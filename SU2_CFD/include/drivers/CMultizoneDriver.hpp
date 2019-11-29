/*!
 * \file CMultizoneDriver.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 7.0.0 "Blackbird"
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

#include "CDriver.hpp"

/*!
 * \class CMultizoneDriver
 * \brief Class for driving zone-specific iterations.
 * \author R. Sanchez, O. Burghardt
 * \version 7.0.0 "Blackbird"
 */
class CMultizoneDriver : public CDriver {
protected:

  bool fsi;
  bool cht;

  unsigned long TimeIter;

  unsigned short *nVarZone;
  su2double **init_res,      /*!< \brief Stores the initial residual. */
            **residual,      /*!< \brief Stores the current residual. */
            **residual_rel;  /*!< \brief Stores the residual relative to the initial. */

  su2double flow_criteria,
            flow_criteria_rel,
            structure_criteria,
            structure_criteria_rel;

  bool *prefixed_motion;     /*!< \brief Determines if a fixed motion is imposed in the config file. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CMultizoneDriver(char* confFile,
             unsigned short val_nZone,
             SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CMultizoneDriver(void);

  /*!
   * \brief [Overload] Launch the computation for multizone problems.
   */
  void StartSolver();

  /*!
   * \brief Preprocess the multizone iteration
   */
  void Preprocess(unsigned long TimeIter);

  /*!
   * \brief Use a corrector step to prevent convergence issues.
   */
  void Corrector(unsigned short val_iZone);

  /*!
   * \brief Run a Block Gauss-Seidel iteration in all physical zones.
   */
  void Run_GaussSeidel();

  /*!
   * \brief Run a Block-Jacobi iteration in all physical zones.
   */
  void Run_Jacobi();

  /*!
   * \brief Update the dual-time solution within multiple zones.
   */
  void Update();

  /*!
   * \brief Output the solution in solution file.
   */
  void Output(unsigned long TimeIter);

  /*!
   * \brief Check the convergence at the outer level.
   */
  bool OuterConvergence(unsigned long OuterIter);

  /*!
   * \brief Perform a dynamic mesh deformation, included grid velocity computation and the update of the multigrid structure (multiple zone).
   */
  void DynamicMeshUpdate(unsigned long TimeIter);

  /*!
   * \brief Perform a dynamic mesh deformation, including grid velocity computation and update of the multigrid structure.
   */
  void DynamicMeshUpdate(unsigned short val_iZone, unsigned long TimeIter);

  /*!
   * \brief Routine to provide all the desired physical transfers between the different zones during one iteration.
   * \return Boolean that determines whether the mesh needs to be updated for this particular transfer
   */
  bool Transfer_Data(unsigned short donorZone, unsigned short targetZone);

  bool Monitor(unsigned long TimeIter);

};
