/*!
 * \class CDiscAdjMultizoneDriver.hpp
 * \brief Class for driving adjoint multi-zone problems.
 * \author O. Burghardt, T. Albring, R. Sanchez
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
#include "CMultizoneDriver.hpp"

class CDiscAdjMultizoneDriver : public CMultizoneDriver {
protected:

  bool retape;                      /*!< \brief Boolean whether derivative information for all zones is kept in memory.*/
  unsigned short RecordingState;    /*!< \brief The kind of recording the tape currently holds.*/
  su2double ObjFunc;                /*!< \brief The value of the objective function.*/
  int ObjFunc_Index;                /*!< \brief Index of the value of the objective function.*/
  unsigned short* direct_nInst;     /*!< \brief Total number of instances in the direct problem (per zone). */
  CIteration*** direct_iteration;   /*!< \brief A pointer to the direct iteration.*/
  COutput** direct_output;          /*!< \brief A pointer to the direct output.*/

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CDiscAdjMultizoneDriver(char* confFile,
             unsigned short val_nZone,
             SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjMultizoneDriver(void);

  /*!
   * \brief [Overload] Launch the computation for discrete adjoint multizone problems.
   */
  void StartSolver();

  /*!
   * \brief Run an discrete adjoint update of all solvers within multiple zones.
   */
  void Run();

  /*!
   * \brief Output the solution in solution file.
   */
  void Output(unsigned long TimeIter);

//  /*!
//   * \brief Check the convergence at the outer level.
//   */
//  bool OuterConvergence(unsigned long OuterIter);

  /*!
   * \brief Record one iteration of a flow iteration in within multiple zones.
   * \param[in] kind_recording - Type of recording (either FLOW_CONS_VARS, MESH_COORDS, COMBINED or NONE)
   * \param[in] indicator which part of a solution update will be recorded
   * \param[in] record_zone - zone where solution update will be recorded
   */
  void SetRecording(unsigned short kind_recording, unsigned short tape_type, unsigned short record_zone);

  /*!
   * \brief Run one direct iteration in a zone.
   */
  void DirectIteration(unsigned short iZone, unsigned short kind_recording);

  /*!
   * \brief Set the objective function. It is virtual because it depends on the kind of physics.
   */
  void SetObjFunction(unsigned short kind_recording);

  /*!
   * \brief Initialize the adjoint value of the objective function.
   */
  void SetAdj_ObjFunction();

  /*!
   * \brief Summary of all routines to evaluate the adjoints in iZone.
   * \param[in] iZone - Zone in which adjoints are evaluated depending on their (preceding) seeding.
   */
  void ComputeAdjoints(unsigned short iZone);

  /*!
   * \brief Saves the current solution (adjoint) values to Solution_Old.
   */
  void Set_OldAdjoints(unsigned short iZone);

  /*!
   * \brief Sets the current iterated solution (adjoint) values to zero.
   */
  void SetIter_Zero(void);

  /*!
   * \brief Adds the current solution (adjoint) values to the iterated solution.
   */
  void Add_IterAdjoints(void);

  /*!
   * \brief Set the current solution (adjoint) values to the current iterated ones.
   */
  void SetAdjoints_Iter(void);

  /*!
   * \brief Computing the RMS residual on driver level (since we iterate zone-wise).
   */
  void SetResidual_RMS(void);
};
