/*!
 * \file COneShotSinglezoneDriver.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
 * \author E. C. Bunschoten
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
#include "CSinglezoneDriver.hpp"
#include "CDiscAdjSinglezoneDriver.hpp"
/*!
 * \class COneShotSinglezoneDriver
 * \ingroup DiscAdj
 * \brief Class for driving single-zone adjoint solvers.
 * \author E. C. Bunschoten
 * \version 8.0.0 "Harrier"
 */
class COneShotSinglezoneDriver : public CDiscAdjSinglezoneDriver {
protected:
    unsigned long n_Iter;                  /*!< \brief The number of adjoint iterations that are run on the fixed-point solver.*/
  RECORDING RecordingState;                     /*!< \brief The kind of recording the tape currently holds.*/
  RECORDING MainVariables;                      /*!< \brief The kind of recording linked to the main variables of the problem.*/
  RECORDING SecondaryVariables;                 /*!< \brief The kind of recording linked to the secondary variables of the problem.*/
  int MainSolver;                               /*!< \brief Index of the main adjoint solver. */
  su2double ObjFunc;                            /*!< \brief The value of the objective function.*/
  CIteration* direct_iteration;                 /*!< \brief A pointer to the direct iteration.*/
  CIteration *iteration;                        /*!< \brief Container vector with all the iteration methods. */
  CIteration *mesh_iteration;

  CConfig *config;                              /*!< \brief Definition of the particular problem. */
  CIntegration **integration;                   /*!< \brief Container vector with all the integration methods. */
  CGeometry *geometry;                          /*!< \brief Geometrical definition of the problem. */
  CSolver **solver;                             /*!< \brief Container vector with all the solutions. */
  COutput *direct_output;
  CNumerics ***numerics;                        /*!< \brief Container vector with all the numerics. */

  su2double** Gradient;
    void ComputeDesignVarGradients();
    void DeformMesh();
    void UpdateDesignVars();
public: 

    /*!
    * \brief Constructor of the class.
    * \param[in] confFile - Configuration file name.
    * \param[in] val_nZone - Total number of zones.
    * \param[in] val_nDim - Total number of dimensions.
    * \param[in] MPICommunicator - MPI communicator for SU2.
    */
    COneShotSinglezoneDriver(char* confFile,
                unsigned short val_nZone,
                SU2_Comm MPICommunicator);
    
    /*!
    * \brief Destructor of the class.
    */
    ~COneShotSinglezoneDriver(void) override;

    void Run(void) override;

    
};