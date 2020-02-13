/*!
 * \file CDiscAdjFEASolver.hpp
 * \brief Headers of the CDiscAdjFEASolver class
 * \author R. Sanchez
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

#include "CSolver.hpp"
#include "../variables/CDiscAdjFEABoundVariable.hpp"

/*!
 * \class CDiscAdjFEASolver
 * \brief Main class for defining the discrete adjoint solver for FE structural problems.
 * \ingroup Discrete_Adjoint
 * \author R. Sanchez
 */
class CDiscAdjFEASolver final : public CSolver {
private:
  unsigned short KindDirect_Solver;
  CSolver *direct_solver;
  su2double *Sens_E,            /*!< \brief Young modulus sensitivity coefficient for each boundary. */
  *Sens_Nu,                     /*!< \brief Poisson's ratio sensitivity coefficient for each boundary. */
  *Sens_nL,                     /*!< \brief Normal pressure sensitivity coefficient for each boundary. */
  **CSensitivity;               /*!< \brief Shape sensitivity coefficient for each boundary and vertex. */

  su2double *Solution_Vel,      /*!< \brief Velocity componenent of the solution. */
  *Solution_Accel;              /*!< \brief Acceleration componenent of the solution. */

  su2double *SolRest;            /*!< \brief Auxiliary vector to restart the solution */

  su2double ObjFunc_Value;      /*!< \brief Value of the objective function. */
  su2double *normalLoads;       /*!< \brief Values of the normal loads for each marker iMarker_nL. */
  unsigned long nMarker_nL;     /*!< \brief Total number of markers that have a normal load applied. */

  unsigned short nMPROP;        /*!< \brief Number of material properties */

  su2double *E_i,               /*!< \brief Values of the Young's Modulus. */
            *Nu_i,              /*!< \brief Values of the Poisson's ratio. */
            *Rho_i,             /*!< \brief Values of the density (for inertial effects). */
            *Rho_DL_i;          /*!< \brief Values of the density (for volume loading). */
  int       *AD_Idx_E_i,        /*!< \brief Derivative index of the Young's Modulus. */
            *AD_Idx_Nu_i,       /*!< \brief Derivative index of the Poisson's ratio. */
            *AD_Idx_Rho_i,      /*!< \brief Derivative index of the density (for inertial effects). */
            *AD_Idx_Rho_DL_i;   /*!< \brief Derivative index of the density (for volume loading). */

  su2double *Local_Sens_E,        /*!< \brief Local sensitivity of the Young's modulus. */
            *Global_Sens_E,       /*!< \brief Global sensitivity of the Young's modulus. */
            *Total_Sens_E;        /*!< \brief Total sensitivity of the Young's modulus (time domain). */
  su2double *Local_Sens_Nu,       /*!< \brief Local sensitivity of the Poisson ratio. */
            *Global_Sens_Nu,      /*!< \brief Global sensitivity of the Poisson ratio. */
            *Total_Sens_Nu;       /*!< \brief Total sensitivity of the Poisson ratio (time domain). */
  su2double *Local_Sens_Rho,      /*!< \brief Local sensitivity of the density. */
            *Global_Sens_Rho,     /*!< \brief Global sensitivity of the density. */
            *Total_Sens_Rho;      /*!< \brief Total sensitivity of the density (time domain). */
  su2double *Local_Sens_Rho_DL,   /*!< \brief Local sensitivity of the volume load. */
            *Global_Sens_Rho_DL,  /*!< \brief Global sensitivity of the volume load. */
            *Total_Sens_Rho_DL;   /*!< \brief Total sensitivity of the volume load (time domain). */

  bool de_effects;                /*!< \brief Determines if DE effects are considered. */
  unsigned short nEField;         /*!< \brief Number of electric field areas in the problem. */
  su2double *EField;              /*!< \brief Array that stores the electric field as design variables. */
  int       *AD_Idx_EField;       /*!< \brief Derivative index of the electric field as design variables. */
  su2double *Local_Sens_EField,   /*!< \brief Local sensitivity of the Electric Field. */
            *Global_Sens_EField,  /*!< \brief Global sensitivity of the Electric Field. */
            *Total_Sens_EField;   /*!< \brief Total sensitivity of the Electric Field (time domain). */

  bool fea_dv;                /*!< \brief Determines if the design variable we study is a FEA parameter. */
  unsigned short nDV;         /*!< \brief Number of design variables in the problem. */
  su2double *DV_Val;          /*!< \brief Values of the design variables. */
  int       *AD_Idx_DV_Val;   /*!< \brief Derivative index of the design variables. */
  su2double *Local_Sens_DV,   /*!< \brief Local sensitivity of the design variables. */
            *Global_Sens_DV,  /*!< \brief Global sensitivity of the design variables. */
            *Total_Sens_DV;   /*!< \brief Total sensitivity of the design variables (time domain). */

  CDiscAdjFEABoundVariable* nodes = nullptr;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

public:

  /*!
   * \brief Constructor of the class.
   */
  CDiscAdjFEASolver(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  CDiscAdjFEASolver(CGeometry *geometry, CConfig *config);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - Initialize the discrete adjoint solver with the corresponding direct solver.
   * \param[in] Kind_Solver - The kind of direct solver.
   */
  CDiscAdjFEASolver(CGeometry *geometry, CConfig *config, CSolver* solver, unsigned short Kind_Solver, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjFEASolver(void);

  /*!
   * \brief Performs the preprocessing of the adjoint AD-based solver.
   *        Registers all necessary variables on the tape. Called while tape is active.
   * \param[in] geometry_container - The geometry container holding all grid levels.
   * \param[in] config_container - The particular config.
   */
  void RegisterSolution(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Performs the preprocessing of the adjoint AD-based solver.
   *        Registers all necessary variables that are output variables on the tape.
   *        Called while tape is active.
   * \param[in] geometry_container - The geometry container holding all grid levels.
   * \param[in] config_container - The particular config.
   */
  void RegisterOutput(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Sets the adjoint values of the output of the flow (+turb.) iteration
   *         before evaluation of the tape.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] config - The particular config.
   */
  void SetAdjoint_Output(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Sets the adjoint values of the input variables of the flow (+turb.) iteration
   *        after tape has been evaluated.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] config - The particular config.
   */
  void ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Sets the adjoint values of the structural variables due to cross term contributions
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  void ExtractAdjoint_CrossTerm(CGeometry *geometry,  CConfig *config) override;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  void ExtractAdjoint_CrossTerm_Geometry(CGeometry *geometry,  CConfig *config) override;

  /*!
   * \brief Register the objective function as output.
   * \param[in] geometry - The geometrical definition of the problem.
   */
  void RegisterObj_Func(CConfig *config) override;

  /*!
   * \brief Set the surface sensitivity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSurface_Sensitivity(CGeometry *geometry, CConfig* config) override;

  /*!
   * \brief Extract and set the geometrical sensitivity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - The solver container holding all terms of the solution.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSensitivity(CGeometry *geometry, CSolver **solver, CConfig *config) override;

  /*!
   * \brief Set the objective function.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetAdj_ObjFunc(CGeometry *geometry, CConfig* config) override;

  /*!
   * \brief Provide the total Young's modulus sensitivity
   * \return Value of the total Young's modulus sensitivity
   *         (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Sens_E(unsigned short iVal) const override { return Total_Sens_E[iVal]; }

  /*!
   * \brief Set the total Poisson's ratio sensitivity.
   * \return Value of the Poisson's ratio sensitivity
   */
  inline su2double GetTotal_Sens_Nu(unsigned short iVal) const override { return Total_Sens_Nu[iVal]; }

  /*!
   * \brief Get the total sensitivity for the structural density
   * \return Value of the structural density sensitivity
   */
  inline su2double GetTotal_Sens_Rho(unsigned short iVal) const override { return Total_Sens_Rho[iVal]; }

  /*!
   * \brief Get the total sensitivity for the structural weight
   * \return Value of the structural weight sensitivity
   */
  inline su2double GetTotal_Sens_Rho_DL(unsigned short iVal) const override { return Total_Sens_Rho_DL[iVal]; }

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Electric Field in the region iEField (time averaged)
   */
  inline su2double GetTotal_Sens_EField(unsigned short iEField) const override { return Total_Sens_EField[iEField]; }

  /*!
   * \brief A virtual member.
   * \return Value of the total sensitivity coefficient for the FEA DV in the region iDVFEA (time averaged)
   */
  inline su2double GetTotal_Sens_DVFEA(unsigned short iDVFEA) const override { return Total_Sens_DV[iDVFEA]; }

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Young Modulus E
   */
  inline su2double GetGlobal_Sens_E(unsigned short iVal) const override { return Global_Sens_E[iVal]; }

  /*!
   * \brief A virtual member.
   * \return Value of the Mach sensitivity for the Poisson's ratio Nu
   */
  inline su2double GetGlobal_Sens_Nu(unsigned short iVal) const override { return Global_Sens_Nu[iVal]; }

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Electric Field in the region iEField
   */
  inline su2double GetGlobal_Sens_EField(unsigned short iEField) const override { return Global_Sens_EField[iEField]; }

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the FEA DV in the region iDVFEA
   */
  inline su2double GetGlobal_Sens_DVFEA(unsigned short iDVFEA) const override { return Global_Sens_DV[iDVFEA]; }

  /*!
   * \brief Get the total sensitivity for the structural density
   * \return Value of the structural density sensitivity
   */
  inline su2double GetGlobal_Sens_Rho(unsigned short iVal) const override { return Global_Sens_Rho[iVal]; }

  /*!
   * \brief Get the total sensitivity for the structural weight
   * \return Value of the structural weight sensitivity
   */
  inline su2double GetGlobal_Sens_Rho_DL(unsigned short iVal) const override { return Global_Sens_Rho_DL[iVal]; }

  /*!
   * \brief Get the value of the Young modulus from the adjoint solver
   * \return Value of the Young modulus from the adjoint solver
   */
  inline su2double GetVal_Young(unsigned short iVal) const override { return E_i[iVal]; }

  /*!
   * \brief Get the value of the Poisson's ratio from the adjoint solver
   * \return Value of the Poisson's ratio from the adjoint solver
   */
  inline su2double GetVal_Poisson(unsigned short iVal) const override { return Nu_i[iVal]; }

  /*!
   * \brief Get the value of the density from the adjoint solver, for inertial effects
   * \return Value of the density from the adjoint solver
   */
  inline su2double GetVal_Rho(unsigned short iVal) const override { return Rho_i[iVal]; }

  /*!
   * \brief Get the value of the density from the adjoint solver, for dead loads
   * \return Value of the density for dead loads, from the adjoint solver
   */
  inline su2double GetVal_Rho_DL(unsigned short iVal) const override { return Rho_DL_i[iVal]; }

  /*!
   * \brief Get the number of variables for the Electric Field from the adjoint solver
   * \return Number of electric field variables from the adjoint solver
   */
  inline unsigned short GetnEField(void) const override { return nEField; }

  /*!
   * \brief Read the design variables for the adjoint solver
   */
  void ReadDV(CConfig *config) override;

  /*!
   * \brief Get the number of design variables from the adjoint solver,
   * \return Number of design variables from the adjoint solver
   */
  inline unsigned short GetnDVFEA(void) const override { return nDV; }

  /*!
   * \brief Get the value of the Electric Field from the adjoint solver
   * \return Pointer to the values of the Electric Field
   */
  inline su2double GetVal_EField(unsigned short iVal) const override { return EField[iVal]; }

  /*!
   * \brief Get the value of the design variables from the adjoint solver
   * \return Pointer to the values of the design variables
   */
  inline su2double GetVal_DVFEA(unsigned short iVal) const override { return DV_Val[iVal]; }

  /*!
   * \brief Prepare the solver for a new recording.
   * \param[in] kind_recording - Kind of AD recording.
   */
  void SetRecording(CGeometry *geometry, CConfig *config) override ;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reset - If true reset variables to their initial values.
   */
  void RegisterVariables(CGeometry *geometry,
                         CConfig *config,
                         bool reset = false) override;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Update the dual-time derivatives.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry,
                    CSolver **solver_container,
                    CConfig *config,
                    unsigned short iMesh,
                    unsigned short iRKStep,
                    unsigned short RunTime_EqSystem,
                    bool Output) override;

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry,
                   CSolver ***solver,
                   CConfig *config,
                   int val_iter,
                   bool val_update_geo) override;

  /*!
   * \brief Compute the multizone residual.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual_Multizone(CGeometry *geometry, CConfig *config) override;

};