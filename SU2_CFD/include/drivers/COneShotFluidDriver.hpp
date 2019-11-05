/*!
 * \file COneShotFluidOutput.hpp
 * \brief Headers of the main subroutines for driving one-shot problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
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
#include "CDiscAdjSinglezoneDriver.hpp"

/*!
 * \class COneShotFluidDriver
 * \brief Class that extends the DiscAdjFluidDriver class with methods for one-shot optimization
 * \author L. Kusch,  B. Munguía
 * \version 6.0.1 "Falcon"
 */
class COneShotFluidDriver : public CDiscAdjSinglezoneDriver {

protected:
  unsigned short RecordingState; /*!< \brief The kind of recording the tape currently holds.*/

  unsigned short nDV_Total; /*!< \brief Total number of design variables used in optimization.*/
  unsigned short nConstr; /*!< \brief Total number of constraints used in optimization.*/

  su2double* Gradient; /*!< \brief Vector to store gradient obtained from projection .*/
  su2double* Gradient_Old; /*!< \brief Vector to store gradient obtained from projection (old value) .*/
  su2double* ShiftedLagrangianGradient; /*!< \brief Saved gradient N_u of the shifted Lagrangian .*/
  su2double* ShiftedLagrangianGradient_Old; /*!< \brief Saved gradient N_u of the shifted Lagrangian (old value) .*/
  su2double* DesignVarUpdate; /*!< \brief Update of the design variable Delta u = u_{k+1} - u_k .*/
  su2double* SearchDirection; /*!< \brief Search direction for optimization.*/
  su2double** BFGS_Inv; /*!< \brief Inverse matrix for BFGS update.*/
  su2double* DesignVariable; /*!< \brief Current design variable value.*/
  su2double* AugmentedLagrangianGradient; /*!< \brief Gradient of doubly augmented Lagrangian.*/
  su2double* AugmentedLagrangianGradient_Old; /*!< \brief Gradient of doubly augmented Lagrangian (old value).*/
  su2double Lagrangian, Lagrangian_Old, Lagrangian_p; /*!< \brief Value of doubly augmented Lagrangian.*/
  su2double GradDotDir; /*!< \brief Gradient dotted with search direction at first Armijo search step (stepsize = 1.0).*/
  
  su2double lb, ub; /*!< \brief Lower and upper bounds of design variables.*/
  su2double epsilon; /*!< \brief Estimator for the active set.*/
  su2double cwolfeone; /*!< \brief First Wolfe line search parameter.*/
  bool* activeset; /*!< \brief Flag for indices belonging to the active set (lower and upper design bounds are reached).*/

  su2double* ConstrFunc; /*!< \brief Constraint function values.*/
  su2double* Multiplier; /*!< \brief Lagrange multipliers for constraint functions.*/
  su2double* Multiplier_Old; /*!< \brief Old Lagrange multipliers for constraint functions.*/
  su2double* ConstrFunc_Store; /*!< \brief Old constraint function (stored when overwritten).*/
  su2double** BCheck_Inv; /*!< \brief Inverse matrix for multiplier update.*/
  su2double  BCheck_Norm; /*!< \brief Norm of the matrix for multiplier update.*/

  //Limited memory BFGS
  unsigned short nBFGSmax; /*!< \brief Maximum number of stored values for BFGS update.*/
  unsigned short nBFGS; /*!< \brief Counter for limited memory BFGS update.*/
  su2double** ykvec; /*!< \brief Vector yk stored for limited memory BFGS.*/
  su2double** skvec; /*!< \brief Vector sk stored for limited memory BFGS.*/

  su2double BFGS_Init;

public:

  /*!
    * \brief Constructor of the class.
    * \param[in] confFile - Configuration file name.
    * \param[in] val_nZone - Total number of zones.
    */
  COneShotFluidDriver(char* confFile,
                   unsigned short val_nZone,
                   SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~COneShotFluidDriver(void);

  /*!
   * \brief Preprocess the one-shot iteration
   * \param[in] TimeIter - index of the current time-step.
   */
  void Preprocess(unsigned long TimeIter);

  /*!
   * \brief Runs main routines of the one-shot class.
   */
  void Run();

  /*!
   * \brief Runs an optimization using the one-shot method.
   */
  void RunOneShot();

  /*!
   * \brief Executes one primal and dual iteration (in a piggy-back manner).
   */
  void PrimalDualStep();

  /*!
   * \brief Executes all operations needed to find the "BCheck"-term of the doubly augmented Lagrangian.
   */
  void ComputeGammaTerm();

  /*!
   * \brief Executes all operations needed to find the "alpha"-term of the doubly augmented Lagrangian.
   */
  void ComputeAlphaTerm();

  /*!
   * \brief Executes all operations needed to find the "beta"-term of the doubly augmented Lagrangian.
   */
  void ComputeBetaTerm();

  /*!
   * \brief Record one iteration of a flow iteration in within multiple zones.
   * \param[in] kind_recording - Type of recording (either CONS_VARS, MESH_COORDS, COMBINED or NONE)
   */
  void SetRecording(unsigned short kind_recording) override;

  /*!
   * \brief Projection of the surface sensitivity using algorithmic differentiation (AD) (see also SU2_DOT).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement class of the problem.
   * \param[in] Gradient - Output to store the gradient data.
   */
  void SetProjection_AD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double* Gradient);

  /*!
   * \brief Performs a surface deformation and volumetric deformation (see also SU2_DEF).
   * \param[in] geometry - geometry class.
   * \param[in] config - config class.
   * \param[in] surface_movement - surface movement class
   * \param[in] grid_movement - volumetric movement class
   */
  void SurfaceDeformation(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, CVolumetricMovement *grid_movement);

  /*!
   * \brief Update the inverse of the BFGS approximation.
   * \param[in] config - config class.
   */
  void BFGSUpdate(CConfig *config);

  /*!
   * \brief Update the inverse of the Limited memory BFGS approximation.
   * \param[in] config - config class.
   */
  void LBFGSUpdate(CConfig *config);

  /*!
   * \brief Recursive application for Limited memory BFGS approximation.
   * \param[in] config - config class.
   */
  void LBFGSUpdateRecursive(CConfig *config, unsigned short nCounter);

  /*!
   * \brief Compute the search direction using the approximated inverse, the gradient N_u and an active set projection.
   */
  void ComputeSearchDirection();
  void ComputeNegativeSearchDirection(); //TODO

  /*!
   * \brief Compute the updated design variable with the search direction and the given step size.
   * \param[in] stepsize - factor for the line search.
   */  
  void ComputeDesignVarUpdate(su2double stepsize);

  /*!
   * \brief Check if the search direction is a descent direction.
   */
  bool CheckDescent();

  /*!
   * \brief Check if the first Wolfe descent condition is fulfilled (line search condition).
   */
  bool CheckFirstWolfe();

  /*!
   * \brief Store gradient dotted with search direction for first Armijo search iteration.
   */
  void StoreGradDotDir();

  /*!
   * \brief Perform parabolic backtracking.
   */
  su2double UpdateStepSizeQuadratic(void);

  /*!
   * \brief Perform cubic backtracking.
   */
  su2double UpdateStepSizeCubic(su2double stepsize, su2double stepsize_p);

  /*!
   * \brief Bound step size.
   */
  su2double UpdateStepSizeBound(su2double stepsize, su2double a, su2double b);

  /*!
   * \brief Store values for the updated design variables.
   */
  void UpdateDesignVariable();

  /*!
   * \brief Store the values for the Lagrangian and its gradient.
   */
  void StoreLagrangianInformation();

  /*!
   * \brief Calculate value for normal or doubly augmented Lagrangian.
   * \param[in] augmented - YES if the doubly augmented Lagrangian shall be calculated
   */
  void CalculateLagrangian(bool augmented);

  /*!
   * \brief Store the gradient of the shifted Lagrangian.
   */
  void SetShiftedLagrangianGradient();

  /*!
   * \brief Store the gradient of the augmented Lagrangian.
   */
  void SetAugmentedLagrangianGradient();

  /*!
   * \brief Initialize the adjoint value of the objective function with 0.0.
   */
  void SetAdj_ObjFunction_Zero();

  /*!
   * \brief Initialize the adjoint value of the constraint function with 0.0.
   */
  void SetAdj_ConstrFunction_Zero();

  /*!
   * \brief Find the indices for the epsilon-active set (overestimated).
   */
  void ComputeActiveSet(su2double stepsize);

  /*!
   * \brief Project mesh sensitivities to surface and design variables using AD.
   */
  void ProjectMeshSensitivities();

  /*!
   * \brief Project input value into feasible bounds of design space.
   * \param[in] value - value that is projected
   * \return projected value
   */
  su2double BoundProjection(su2double value);

  /*!
   * \brief Project a given vector value into the epsilon-active set (or inactive set).
   * \param[in] iDV - index of the respective value
   * \param[in] value - corresponding value
   * \param[in] active - is true if the projection is into the active set
   * \return projected value
   */
  su2double ProjectionSet(unsigned short iDV, su2double value, bool active);

  /*!
   * \brief Project a given matrix value into the epsilon-active set (or inactive set).
   * \param[in] iDV - row index of the respective value
   * \param[in] jDV - column index of the respective value
   * \param[in] value - corresponding value
   * \param[in] active - is true if the projection is into the active set
   * \return projected value
   */
  su2double ProjectionPAP(unsigned short iDV, unsigned short jDV, su2double value, bool active);

  /*!
   * \brief Store the old constraint multiplier.
   */
  void StoreMultiplier();

    /*!
   * \brief Load the old constraint multiplier.
   */
  void LoadMultiplier();

  /*!
   * \brief Update the constraint multiplier.
   */
  void UpdateMultiplier(su2double stepsize);

  /*!
   * \brief Check the sign of the constraint multiplier.
   */
  void CheckMultiplier();

  /*!
   * \brief Set the constraint functions.
   */
  void SetConstrFunction();

  /*!
   * \brief Initialize the adjoint value of the constraint functions.
   */
  void SetAdj_ConstrFunction(su2double *seeding);

  /*!
   * \brief Store the constraints.
   */
  void StoreConstrFunction();

  /*!
   * \brief Compute the inverse preconditioner matrix (BCheck^(-1)) for the multiplier update.
   */
  void ComputePreconditioner();
};