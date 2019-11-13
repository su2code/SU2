/*!
 * \file CMultizoneDriver.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
 * \author T. Economon, H. Kline, R. Sanchez
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

#include "CDriver.hpp"

/*!
 * \class CMultizoneDriver
 * \brief Class for driving zone-specific iterations.
 * \author R. Sanchez, O. Burghardt
 * \version 6.0.1 "Falcon"
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
