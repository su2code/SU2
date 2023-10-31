/*!
 * \file CDiscAdjFEASolver.hpp
 * \brief Headers of the CDiscAdjFEASolver class
 * \author R. Sanchez
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

#include "CSolver.hpp"
#include "../variables/CDiscAdjFEABoundVariable.hpp"

/*!
 * \class CDiscAdjFEASolver
 * \ingroup DiscAdj
 * \brief Main class for defining the discrete adjoint solver for FE structural problems.
 * \author R. Sanchez
 */
class CDiscAdjFEASolver final : public CSolver {
private:
  static constexpr size_t MAXNVAR = 9;  /*!< \brief Max number of variables, for static arrays. */
  unsigned short KindDirect_Solver = 0;
  CSolver *direct_solver = nullptr;

  /*!
   * \brief A type to manage sensitivities of design variables.
   */
  struct SensData {
    unsigned short size = 0;
    su2double* val = nullptr;         /*!< \brief Value of the variable. */
    su2double* LocalSens = nullptr;   /*!< \brief Local sensitivity (domain). */
    su2double* GlobalSens = nullptr;  /*!< \brief Global sensitivity (mpi). */
    su2double* OldSens = nullptr;     /*!< \brief Previous global sensitivity, used to update the total. */
    su2double* TotalSens = nullptr;   /*!< \brief Total sensitivity (integrated over time). */

    su2double& operator[] (unsigned short i) { return val[i]; }
    const su2double& operator[] (unsigned short i) const { return val[i]; }

    void resize(unsigned short n) {
      clear();
      size = n;
      val = new su2double[n]();
      LocalSens = new su2double[n]();
      GlobalSens = new su2double[n]();
      OldSens = new su2double[n]();
      TotalSens = new su2double[n]();
    }

    void clear() {
      size = 0;
      delete [] val;
      delete [] LocalSens;
      delete [] GlobalSens;
      delete [] OldSens;
      delete [] TotalSens;
    }

    void Register() {
      for (auto i = 0u; i < size; ++i) AD::RegisterInput(val[i]);
    }

    void GetDerivative() {
      for (auto i = 0u; i < size; ++i) LocalSens[i] = SU2_TYPE::GetDerivative(val[i]);

      SU2_MPI::Allreduce(LocalSens, GlobalSens, size, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

      for (auto i = 0u; i < size; ++i) {
        /*--- Update the total by subtracting the old and adding the new value.
         * Then update the old value for the next call to this function. ---*/
        TotalSens[i] += GlobalSens[i] - OldSens[i];
        OldSens[i] = GlobalSens[i];
      }
    }

    void Store() {
      /*--- Clears the old values such that on the next time step the total is
       * incremented instead of updated. ---*/
      for (auto i = 0u; i < size; ++i) OldSens[i] = 0.0;
    }

    ~SensData() { clear(); }
  };

  unsigned short nMPROP = 0;  /*!< \brief Number of material properties */
  SensData E;                 /*!< \brief Values of the Young's Modulus. */
  SensData Nu;                /*!< \brief Values of the Poisson's ratio. */
  SensData Rho;               /*!< \brief Values of the density (for inertial effects). */
  SensData Rho_DL;            /*!< \brief Values of the density (for volume loading). */

  bool de_effects = false;    /*!< \brief Determines if DE effects are considered. */
  unsigned short nEField = 0; /*!< \brief Number of electric field areas in the problem. */
  SensData EField;            /*!< \brief Array that stores the electric field as design variables. */

  bool fea_dv = false;        /*!< \brief Determines if the design variable we study is a FEA parameter. */
  unsigned short nDV = 0;     /*!< \brief Number of design variables in the problem. */
  SensData DV;                /*!< \brief Values of the design variables. */

  CDiscAdjFEABoundVariable* nodes = nullptr;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

  /*!
   * \brief Read the design variables for the adjoint solver
   */
  void ReadDV(const CConfig *config);

public:

  /*!
   * \brief Constructor of the class.
   */
  CDiscAdjFEASolver() = default;

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
  ~CDiscAdjFEASolver() override;

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
  void ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config, bool CrossTerm) override;

  /*!
   * \brief Extract and set the geometrical sensitivity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSensitivity(CGeometry *geometry, CConfig *config, CSolver*) override;

  /*!
   * \brief Provide the total Young's modulus sensitivity
   * \return Value of the total Young's modulus sensitivity
   *         (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Sens_E(unsigned short iVal) const override { return E.TotalSens[iVal]; }

  /*!
   * \brief Set the total Poisson's ratio sensitivity.
   * \return Value of the Poisson's ratio sensitivity
   */
  inline su2double GetTotal_Sens_Nu(unsigned short iVal) const override { return Nu.TotalSens[iVal]; }

  /*!
   * \brief Get the total sensitivity for the structural density
   * \return Value of the structural density sensitivity
   */
  inline su2double GetTotal_Sens_Rho(unsigned short iVal) const override { return Rho.TotalSens[iVal]; }

  /*!
   * \brief Get the total sensitivity for the structural weight
   * \return Value of the structural weight sensitivity
   */
  inline su2double GetTotal_Sens_Rho_DL(unsigned short iVal) const override { return Rho_DL.TotalSens[iVal]; }

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Electric Field in the region iEField (time averaged)
   */
  inline su2double GetTotal_Sens_EField(unsigned short iEField) const override { return EField.TotalSens[iEField]; }

  /*!
   * \brief A virtual member.
   * \return Value of the total sensitivity coefficient for the FEA DV in the region iDVFEA (time averaged)
   */
  inline su2double GetTotal_Sens_DVFEA(unsigned short iDVFEA) const override { return DV.TotalSens[iDVFEA]; }

  /*!
   * \brief Get the value of the Young modulus from the adjoint solver
   * \return Value of the Young modulus from the adjoint solver
   */
  inline su2double GetVal_Young(unsigned short iVal) const override { return E[iVal]; }

  /*!
   * \brief Get the value of the Poisson's ratio from the adjoint solver
   * \return Value of the Poisson's ratio from the adjoint solver
   */
  inline su2double GetVal_Poisson(unsigned short iVal) const override { return Nu[iVal]; }

  /*!
   * \brief Get the value of the density from the adjoint solver, for inertial effects
   * \return Value of the density from the adjoint solver
   */
  inline su2double GetVal_Rho(unsigned short iVal) const override { return Rho[iVal]; }

  /*!
   * \brief Get the value of the density from the adjoint solver, for dead loads
   * \return Value of the density for dead loads, from the adjoint solver
   */
  inline su2double GetVal_Rho_DL(unsigned short iVal) const override { return Rho_DL[iVal]; }

  /*!
   * \brief Get the number of variables for the Electric Field from the adjoint solver
   * \return Number of electric field variables from the adjoint solver
   */
  inline unsigned short GetnEField(void) const override { return nEField; }

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
  inline su2double GetVal_DVFEA(unsigned short iVal) const override { return DV[iVal]; }

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

};
