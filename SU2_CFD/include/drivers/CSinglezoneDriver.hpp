/*!
 * \file CSinglezoneDriver.hpp
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
 * \class CSinglezoneDriver
 * \brief Class for driving single-zone solvers.
 * \author R. Sanchez
 * \version 7.0.0 "Blackbird"
 */
class CSinglezoneDriver : public CDriver {
protected:

  unsigned long TimeIter;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CSinglezoneDriver(char* confFile,
             unsigned short val_nZone,
             SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CSinglezoneDriver(void);

  /*!
   * \brief [Overload] Launch the computation for single-zone problems.
   */
  void StartSolver();

  /*!
   * \brief Preprocess the single-zone iteration
   */
  virtual void Preprocess(unsigned long TimeIter);

  /*!
   * \brief Run the iteration for ZONE_0.
   */
  virtual void Run();

  /*!
   * \brief Postprocess the iteration for ZONE_0.
   */
  virtual void Postprocess();

  /*!
   * \brief Update the dual-time solution within multiple zones.
   */
  void Update();

  /*!
   * \brief Output the solution in solution file.
   */
  void Output(unsigned long TimeIter);

  /*!
   * \brief Perform a dynamic mesh deformation, included grid velocity computation and the update of the multigrid structure.
   */
  void DynamicMeshUpdate(unsigned long TimeIter);

  /*!
   * \brief Monitor
   * \param ExtIter
   */
  virtual bool Monitor(unsigned long TimeIter);

  /*!
     * \brief  Returns wheter all specified windowed-time-averaged ouputs have been converged
     * \return Boolean indicating whether the problem is converged.
     */
  inline virtual bool GetTimeConvergence() const{
    return output_container[ZONE_0]->GetTimeConvergence();
  }


  /*!
   * \brief Runtime_Parsing
   */
  virtual void Runtime_Options();

};
