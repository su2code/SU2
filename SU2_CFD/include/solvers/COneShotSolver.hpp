/*!
 * \file COneShotSolver.hpp
 * \brief Main header for defining the OneShot solver.
 * \author L. Kusch, T.Dick
 * \version 7.0.6 "Blackbird"
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

#include "CDiscAdjSolver.hpp"

class COneShotSolver : public CDiscAdjSolver {
private:
  su2double theta, rho;
  unsigned short nConstr;
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
   * \brief Prepare the solver for a new recording (without setting solution to initial solution).
   * \param[in] geometry - geometry class object
   * \param[in] config - config class object
   */
  void SetRecording(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Store the current solution in Solution_Store.
   * (This is done to store the solution of the current iterate before line search)
   */
  void StoreSolution();

  /*!
   * \brief Load the current solution from Solution_Store.
   */
  void LoadSolution();

  /*!
   * \brief Load the current solution from Solution_Store + stepsize*Update.
   */
  void LoadSolutionStep(su2double stepsize);

  /*!
   * \brief Store current mesh coordinates and normals.
   * (This is e.g. done before line search)
   * \param[in] config - config class object
   * \param[in] geometry - geometry class object
   */
  void StoreMeshPoints(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Load mesh coordinates and normals.
   * \param[in] config - config class object
   * \param[in] geometry - geometry class object
   */
  void LoadMeshPoints(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Store the geometry sensitivities (Sensitivity) in Sensitivity_ShiftedLagrangian
   * \param[in] geometry - geometry class object
   */
  void SaveSensitivity(CGeometry* geometry);

  /*!
   * \brief Calculate either the solver part of the augmented or the shifted Lagrangian
   * (without objective and constraint functions)
   * \param[in] config - config class object
   * \param[in] augmented - true if the augmented Lagrangian is considered
   * \result value of the Lagrangian part
   */
  su2double CalculateLagrangianPart(CConfig* config, bool augmented);

  /*!
   * \brief Store the current solution in Solution_Save.
   * (This is done to store the solution before calculating the alpha and beta terms)
   */
  void StoreSaveSolution();

  /*!
   * \brief Load the current solution from Solution_Save.
   */
  void LoadSaveSolution();

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
  void SetGeometrySensitivityLagrangian(CGeometry *geometry);

  /*!
   * \brief Shift Solution_Former to Solution_Store.
   * (This is needed if the solution at two former states is stored)
   */
  void ShiftFormerSolution();

  /*!
   * \brief Shift Solution_Store to Solution_Former.
   * (This is needed if Solution is only updated by a factor*Delta)
   */
  void ShiftStoreSolution();

  /*!
   * \brief Store Solution_Delta in Solution_Delta_Store.
   * (This is needed if Solution is only updated by a factor*Delta)
   */
  void StoreSolutionDelta();

  /*!
   * \brief Store the current solution in Solution_Former.
   * (This additional storage is needed if the solution at two former states is stored)
   */
  void StoreFormerSolution();

  /*!
   * \brief Calculate estimates for alpha, beta, and gamma of the doubly augmented Lagrangian
   */
  void CalculateAlphaBetaGamma(CConfig *config);

  /*!
   * \brief Store estimates for alpha, beta, and gamma of the doubly augmented Lagrangian
   */
  void SetAlphaBetaGamma(CConfig *config, su2double val_bcheck_norm);

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
   * \brief Set the sensitivity of the doubly augmented Lagrangian to zero.
   * \param[in] geometry - geometry class element
   */
  void ResetSensitivityLagrangian(CGeometry *geometry);

  /*!
   * \brief Update the sensitivity of the doubly augmented Lagrangian with a factor.
   * \param[in] geometry - geometry class element
   * \param[in] factor - multiplier for the update
   */
  void UpdateSensitivityLagrangian(CGeometry *geometry, su2double factor);

  /*!
   * \brief Adds the difference in the adjoint solution to the solution (with a finite difference step size).
   * (This is done to calculate the finite difference update for the beta term)
   * \param[in] config - config class element
   */
  void UpdateStateVariable(CConfig *config);

  /*!
   * \brief Set the sensitivity to the finite difference to approximate N_yx for the beta term.
   * \param[in] geometry - geometry class element
   * \param[in] config - config class element
   */
  void SetFiniteDifferenceSens(CGeometry *geometry, CConfig *config);

  void SetConstrDerivative(unsigned short iConstr);

  su2double MultiplyConstrDerivative(unsigned short iConstr, unsigned short jConstr);

  void SetObjFunc_Value(su2double val_ObjFunc);

  su2double GetObjFunc_Value(void);

  void SetConFunc_Value(su2double val_ConFunc);

  void AddConFunc_Value(su2double val_ConFunc);

  su2double GetConFunc_Value(void);

};
