/*!
 * \file error_estimation_structure.hpp
 * \brief Headers of the main subroutines for creating the error estimation structure.
 *        The subroutines and functions are in the <i>error_estimation_structure.cpp</i> file.
 * \author B. Munguía
 * \version 6.1.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

#include "../../Common/include/adt_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/fem_geometry_structure.hpp"
#include "../../SU2_CFD/include/iteration_structure.hpp"
#include "../../SU2_CFD/include/integration_structure.hpp"
#include "../../SU2_CFD/include/output_structure.hpp"
#include "../../SU2_CFD/include/numerics_structure.hpp"
#include "../../SU2_CFD/include/solver_structure.hpp"

using namespace std;

/*!
 * \class CErrorEstimationDriver
 * \brief Class to drive error estimation.
 * \version 6.1.0 "Falcon"
 */
class CErrorEstimationDriver {
private:
  int rank,                                     /*!< \brief MPI Rank. */
      size;                                     /*!< \brief MPI Size. */

  COutput *output;                              /*!< \brief Pointer to the COutput class. */
  CGeometry ****geometry_container;             /*!< \brief Geometry for which there is a solution */
  CSolver *****solver_container;                /*!< \brief Solution for which error is being estimated */
  CConfig **config_container;                   /*!< \brief Definition of the problem. */
  CIteration ***iteration_container;            /*!< \brief Container vector with all the iteration methods. */
  CIntegration ****integration_container;       /*!< \brief Container vector with all the integration methods. */
  CNumerics ******numerics_container;           /*!< \brief Description of the numerical method (the way in which the equations are solved). */
  CConfig *driver_config;                       /*!< \brief Definition of the driver configuration. */
  CSurfaceMovement **surface_movement;          /*!< \brief Surface movement classes of the problem. */
  CVolumetricMovement ***grid_movement;         /*!< \brief Volume grid movement classes of the problem. */
  CFreeFormDefBox*** FFDBox;                    /*!< \brief FFD FFDBoxes of the problem. */
  char* config_file_name;                       /*!< \brief Configuration file name of the problem.*/

  bool fsi,                                     /*!< \brief FSI simulation flag.*/
       fem_solver;                              /*!< \brief FEM fluid solver simulation flag */

  unsigned short iZone,                         /*!< \brief Iterator on zones.*/
                 iSol,                          /*!< \brief Iterator on solutions.*/
                 nZone,                         /*!< \brief Total number of zones in the problem. */
                 nDim,                          /*!< \brief Number of dimensions.*/
                 iInst,                         /*!< \brief Iterator on instance levels.*/
                 *nInst;                        /*!< \brief Total number of instances in the problem (per zone). */
  
  unsigned long DOFsPerPoint;                   /*!< \brief Number of unknowns at each vertex, i.e., number of equations solved. */

  unsigned long ExtIter;                        /*!< \brief External iteration.*/

  unsigned short RecordingState;                /*!< \brief The kind of recording the tape currently holds.*/
  su2double ObjFunc;                            /*!< \brief The value of the objective function.*/
  CIteration** direct_iteration;                /*!< \brief A pointer to the direct iteration.*/

  vector<su2double> epsilon_fine;               /*!< \brief A vector to store the adaptation parameter at each fine mesh node. */

public:

  /*! 
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] val_nDim - Number of dimensions.   
   * \param[in] val_periodic - Bool for periodic BCs.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CErrorEstimationDriver(char* confFile,
                          unsigned short val_nZone,
                          unsigned short val_nDim,
                          bool val_periodic,
                          SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CErrorEstimationDriver(void);

  /*!
   * \brief Read in the config and mesh files.
   */
  void Input_Preprocessing(CConfig **config_container, CGeometry ****geometry_container, bool val_periodic);

  /*!
   * \brief Construction of the edge-based data structure and the multigrid structure.
   */
  void Geometrical_Preprocessing(CConfig **config_container, CGeometry ****geometry_container);
  
  /*!
   * \brief Do the geometrical preprocessing for the DG FEM solver.
   */
  void Geometrical_Preprocessing_DGFEM(CConfig **config_container, CGeometry ****geometry_container);

  /*!
   * \brief Construction of iteration classes.
   */
  void Iteration_Preprocessing(void);

  /*!
   * \brief Construction of solution classes.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Solver_Preprocessing(CSolver ****solver_container, CGeometry ***geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Construction of integration classes.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Integration_Preprocessing(CIntegration ***integration_container, CGeometry ***geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Definition and allocation of all numerics classes.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Numerics_Preprocessing(CNumerics *****numerics_container, CSolver ****solver_container, CGeometry ***geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Restart of the solvers from the restart files.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Solver_Restart(CSolver ****solver_container, CGeometry ***geometry, CConfig *config, bool update_geo, unsigned short val_iInst);

  /*!
   * \brief Perform all steps to compute the metric.
   */
  void ComputeMetric(void);

  /*!
   * \brief Perform inner product of adjoint gradients and flux Hessian to compute the adaptation parameter.
   */
  void SumWeightedHessian(CSolver* solver_flow, CSolver* solver_adj);

  /*!
   * \brief Run an iteration of the adjoint solver.
   */
  void Run(void);

  /*!
   * \brief Specify the kind of recording for the adjoint solver.
   */
  void SetRecording(unsigned short kind_recording);

  /*!
   * \brief Set the adjoint value of the objective function.
   */
  void SetAdj_ObjFunction(void);

  /*!
   * \brief Set the objective function.
   */
  void SetObjFunction(void);

  /*!
   * \brief Run an iteration of the flow solver.
   */
  void DirectRun(void);

  /*!
   * \brief Output the solution in solution file.
   */
  void Output(void);

  /*!
   * \brief Deallocation routine
   */
  void Postprocessing(void);

  /*!
   * \brief Deallocation of solution classes.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Solver_Postprocessing(CSolver ****solver_container, CGeometry **geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Deallocation of integration classes.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Integration_Postprocessing(CIntegration ***integration_container, CGeometry **geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Deallocation of solver classes.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Numerics_Postprocessing(CNumerics *****numerics_container, CSolver ***solver_container, CGeometry **geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Definition and allocation of all solution classes.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Solver_Deletion(CSolver ****solver_container, CConfig *config, unsigned short val_iInst);

};