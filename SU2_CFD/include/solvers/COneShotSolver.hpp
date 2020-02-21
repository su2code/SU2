/*!
 * \file COneShotSolver.hpp
 * \brief Headers of the COneShotSolver class
 * \author T. Albring
 * \version 7.0.1 "Blackbird"
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

#include "CDiscAdjSolver.hpp"

/*!
 * \class COneShotSolver
 * \brief Main class for defining the one-shot optimization.
 * \ingroup One_Shot
 * \author L. Kusch
 * \version 5.0.0 "Raven"
 */
class COneShotSolver : public CDiscAdjSolver {
private:
  su2double theta, theta_old, rho, rho_old, grad_norm, *lambda;
  unsigned short nConstr, ArmijoIter, nActiveDV;
  su2double ConFunc_Value;
  su2double *** DConsVec;

public:

  /*!
   * \brief Constructor of the class.
   */
  COneShotSolver(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  COneShotSolver(CGeometry *geometry, CConfig *config);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - Initialize the discrete adjoint solver with the corresponding direct solver.
   * \param[in] Kind_Solver - The kind of direct solver.
   */
  COneShotSolver(CGeometry *geometry, CConfig *config, CSolver* solver, unsigned short Kind_Solver, unsigned short iMesh);

  ~COneShotSolver(void);

  /*!
   * \brief Reset indices of inputs.
   * \param[in] geometry - geometry class object
   * \param[in] config - config class object
   */
  void ResetInputs(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Prepare the solver for a new recording (without setting solution to initial solution).
   * \param[in] geometry - geometry class object
   * \param[in] config - config class object
   */
  void SetRecording(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Store the current solution in Solution_Store.
   * (This is done to store the solution of the current iterate before line search)
   */
  void SetStoreSolution();

  /*!
   * \brief Load the current solution from Solution_Store.
   */
  void LoadSolution();

  /*!
   * \brief Store current mesh coordinates and normals.
   * (This is e.g. done before line search)
   * \param[in] config - config class object
   * \param[in] geometry - geometry class object
   */
  void SetMeshPointsOld(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Load mesh coordinates and normals.
   * \param[in] config - config class object
   * \param[in] geometry - geometry class object
   */
  void LoadMeshPointsOld(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Load old mesh coordinates and add perturbation.
   * \param[in] config - config class object
   * \param[in] geometry - geometry class object
   */
  void LoadMeshPointsStep(CConfig *config, CGeometry *geometry, su2double stepsize);

  /*!
   * \brief Store the geometry sensitivities (Sensitivity) in Sensitivity_ShiftedLagrangian
   * \param[in] geometry - geometry class object
   */
  void SetSensitivityShiftedLagrangian(CGeometry* geometry);

  /*!
   * \brief Calculate either the solver part of the augmented or the shifted Lagrangian
   * (without objective and constraint functions)
   * \param[in] config - config class object
   * \result value of the Lagrangian part
   */
  su2double CalculateLagrangian(CConfig* config);

  /*!
   * \brief Store the current solution in Solution_Save.
   * (This is done to store the solution before calculating the alpha and beta terms)
   */
  void SetSaveSolution();

  /*!
   * \brief Store the previous solution in Solution_Save.
   * (This is done to store the solution before calculating the alpha and beta terms)
   */
  void SetOldStoreSolution();

  /*!
   * \brief Load the current solution from Solution_Save.
   */
  void LoadSaveSolution();

  /*!
   * \brief Increment the solution by stepsize*d(bar)y.
   */
  void LoadStepSolution(su2double stepsize);

  /*!
   * \brief Set the geometry sensitivity to the sensitivity of the shifted Lagrangian
   * (This happens for the sensitivity projection)
   * \param[in] geometry - geometry class object
   */
  void SetGeometrySensitivityGradient(CGeometry *geometry);

  /*!
   * \brief Set the geometry sensitivity to the sensitivity of the augmented Lagrangian
   * (This happens for the sensitivity projection)
   * \param[in] geometry - geometry class object
   */
  void SetGeometrySensitivityLagrangian(CGeometry *geometry, unsigned short kind);

  /*!
   * \brief Set Solution_Delta for this time step.
   */
  void SetSolutionDelta(CGeometry *geometry);

  /*!
   * \brief Set Solution_Delta for this time step.
   */
  void SetSaveSolutionDelta(CGeometry *geometry);

  /*!
   * \brief Set Solution_Delta for this time step.
   */
  void SetStoreSolutionDelta();

  /*!
   * \brief Calculate estimates for alpha, beta, and gamma of the doubly augmented Lagrangian
   */
  void CalculateRhoTheta(CConfig *config);

  /*!
   * \brief Store estimates for alpha and beta of the doubly augmented Lagrangian
   */
  void CalculateAlphaBeta(CConfig *config);

  /*!
   * \brief Store estimates for gamma of the doubly augmented Lagrangian
   */
  void CalculateGamma(CConfig *config, su2double val_bcheck_norm, su2double* val_constr_func, su2double* val_gamma);

  /*!
   * \brief Sets the adjoint values of the input variables of the flow (+turb.) iteration
   *        after tape has been evaluated without computing the adjoint residual.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] config - The particular config.
   */
  void ExtractAdjoint_Solution_Clean(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Set the adjoint output to the difference in the solution.
   * (Is done to evaluate the adjoint solution for the alpha term)
   * \param[in] geometry - geometry class element
   * \param[in] config - config class element
   */
  void SetAdjoint_OutputUpdate(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Set the adjoint output to zero.
   * (Is done to evaluate the adjoint solution for the constraint function only)
   * \param[in] geometry - geometry class element
   * \param[in] config - config class element
   */
  void SetAdjoint_OutputZero(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Update the sensitivity of the doubly augmented Lagrangian with a factor.
   * \param[in] geometry - geometry class element
   * \param[in] factor - multiplier for the update
   */
  void SetSensitivityLagrangian(CGeometry *geometry, unsigned short kind);

  /*!
   * \brief Adds the difference in the adjoint solution to the solution (with a finite difference step size).
   * (This is done to calculate the finite difference update for the beta term)
   * \param[in] config - config class element
   */
  void UpdateStateVariable(CConfig *config, su2double fd_step);

  /*!
   * \brief Set the sensitivity to the finite difference to approximate N_yx for the beta term.
   * \param[in] geometry - geometry class element
   * \param[in] config - config class element
   */
  void SetFiniteDifferenceSens(CGeometry *geometry, CConfig *config);

  void SetConstrDerivative(unsigned short iConstr);

  su2double GetConstrDerivative(unsigned short iConstr, unsigned long iPoint, unsigned long iVar) const { return DConsVec[iConstr][iPoint][iVar]; }

  su2double MultiplyConstrDerivative(unsigned short iConstr, unsigned short jConstr);

  void SetConFunc_Value(su2double val_ConFunc) { ConFunc_Value = val_ConFunc; }

  void AddConFunc_Value(su2double val_ConFunc) { ConFunc_Value += val_ConFunc; }

  su2double GetConFunc_Value(void) const { return ConFunc_Value; }

  void SetArmijoIter(unsigned short val_iter) { ArmijoIter = val_iter; }

  unsigned short GetArmijoIter(void) const { return ArmijoIter; }

  void SetnActiveDV(unsigned short val_active) { nActiveDV = val_active; }

  unsigned short GetnActiveDV(void) const { return nActiveDV; }

  su2double GetOneShotRho(void) { return rho; }

  su2double GetOneShotTheta(void) { return theta; }

  void SetShiftedLagGradNorm(su2double val_norm) { grad_norm = val_norm; }

  su2double GetShiftedLagGradNorm(void) const { return grad_norm; }

  void SetLambdaValue(unsigned short iConstr, su2double val_constr) {lambda[iConstr] = val_constr; }

  su2double GetLambdaValue(unsigned short iConstr) const { return lambda[iConstr]; }

};