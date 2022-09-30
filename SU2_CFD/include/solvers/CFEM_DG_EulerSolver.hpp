/*!
 * \file CFEM_DG_EulerSolver.hpp
 * \brief Headers of the CFEM_DG_EulerSolver class
 * \author E. van der Weide, T. Economon, J. Alonso
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "CFEM_DG_SolverBase.hpp"

/*!
 * \class CFEM_DG_EulerSolver
 * \brief Main class for defining the Euler Discontinuous Galerkin finite element flow solver.
 * \ingroup Euler_Equations
 * \author E. van der Weide, T. Economon, J. Alonso
 * \version 7.4.0 "Blackbird"
 */
class CFEM_DG_EulerSolver : public CFEM_DG_SolverBase {
protected:
  static constexpr size_t MAXNDIM = 3; /*!< \brief Max number of space dimensions, used in some static arrays. */
  static constexpr size_t MAXNVAR = 5; /*!< \brief Max number of variables, for static arrays. */

  su2double Gamma;           /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One; /*!< \brief Fluids's Gamma - 1.0  . */

  vector<su2double> EntropyVarFreeStream; /*!< \brief Vector, which contains the free stream
                                                      entropy variables. */

  vector<unsigned long> nDOFsPerRank;                    /*!< \brief Number of DOFs per rank in
                                                                     cumulative storage format. */
  vector<vector<unsigned long> > nonZeroEntriesJacobian; /*!< \brief The ID's of the DOFs for the
                                                                     non-zero entries of the Jacobian
                                                                     for the locally owned DOFs. */

  int nGlobalColors;              /*!< \brief Number of global colors for the Jacobian computation. */

  vector<vector<unsigned long> > localDOFsPerColor;   /*!< \brief Double vector, which contains for every
                                                                  color the local DOFs. */
  vector<vector<int> > colorToIndEntriesJacobian;     /*!< \brief Double vector, which contains for every
                                                                  local DOF the mapping from the color to the
                                                                  entry in the Jacobian. A -1 indicates that
                                                                  the color does not contribute to the Jacobian
                                                                  of the DOF. */
public:

  /*!
   * \brief Constructor of the class.
   */
  CFEM_DG_EulerSolver(void);

  /*!
   * \overload
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nDim - Dimension of the problem (2D or 3D).
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  CFEM_DG_EulerSolver(CConfig *config, unsigned short val_nDim, unsigned short iMesh);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  CFEM_DG_EulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_DG_EulerSolver(void) override;

  /*!
   * \brief Set the fluid solver nondimensionalization.
   * \param[in] config      - Definition of the particular problem.
   * \param[in] iMesh       - Index of the mesh in multigrid computations.
   * \param[in] writeOutput - Whether or not output must be written.
   */
  void SetNondimensionalization(CConfig        *config,
                                unsigned short iMesh,
                                const bool     writeOutput);

  /*!
   * \brief Set the initial condition for the Euler Equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void SetInitialCondition(CGeometry **geometry,
                           CSolver ***solver_container,
                           CConfig *config,
                           unsigned long TimeIter) final;

  /*!
   * \brief Set the working solution of the first time level to the current
   *        solution. Used for Runge-Kutta type schemes.
   */
  void Set_OldSolution() final;

  /*!
   * \brief Set the new solution to the current solution for classical RK.
   */
  void Set_NewSolution() final;

  /*!
   * \brief Function to compute the time step for solving the Euler equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Value of the current iteration.
   */
  void SetTime_Step(CGeometry      *geometry,
                    CSolver        **solver_container,
                    CConfig        *config,
                    unsigned short iMesh,
                    unsigned long  Iteration) override;

  /*!
   * \brief Function, which checks whether or not the time synchronization point is reached
            when explicit time stepping is used.
   * \param[in]     config          - Definition of the particular problem.
   * \param[in]     TimeSync        - The synchronization time.
   * \param[in,out] timeEvolved     - On input the time evolved before the time step,
                                      on output the time evolved after the time step.
   * \param[out]    syncTimeReached - Whether or not the synchronization time is reached.
   */
  void CheckTimeSynchronization(CConfig         *config,
                                const su2double TimeSync,
                                su2double       &timeEvolved,
                                bool            &syncTimeReached) final;

  /*!
   * \brief Function, which processes the list of tasks to be executed by
            the DG solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void ProcessTaskList_DG(CGeometry      *geometry,
                          CSolver        **solver_container,
                          CNumerics      **numerics,
                          CConfig        *config,
                          unsigned short iMesh) final;

  /*!
   * \brief Function, to carry out the space time integration for ADER
            with time accurate local time stepping.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void ADER_SpaceTimeIntegration(CGeometry      *geometry,
                                 CSolver        **solver_container,
                                 CNumerics      **numerics,
                                 CConfig        *config,
                                 unsigned short iMesh,
                                 unsigned short RunTime_EqSystem) final;

  /*!
   * \brief Function, which controls the computation of the spatial Jacobian.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void ComputeSpatialJacobian(CGeometry      *geometry,
                              CSolver        **solver_container,
                              CNumerics      **numerics,
                              CConfig        *config,
                              unsigned short iMesh,
                              unsigned short RunTime_EqSystem) final;

  /*!
   * \brief Function, which determines the values of the tolerances in
            the predictor step of ADER-DG.
   */
  void TolerancesADERPredictorStep(void);

  /*!
   * \brief Function, carries out the predictor step of the ADER-DG
            time integration.
   * \param[in]  config    - Definition of the particular problem.
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   */
  void ADER_DG_PredictorStep(CConfig             *config,
                             const unsigned long elemBeg,
                             const unsigned long elemEnd);

  /*!
   * \brief Function, which interpolates the predictor solution of ADER-DG
            to the time value that corresponds to iTime.
   * \param[in]  config            - Definition of the particular problem.
   * \param[in]  iTime             - Time index of the time integration point for the
                                     integration over the time slab in the corrector
                                     step of ADER-DG.
   * \param[in]  elemBeg           - Begin index of the element range to be computed. This
                                     range is for elements of the same time level.
   * \param[in]  elemEnd           - End index (not included) of the element range to be computed.
   * \param[in]  nAdjElem          - Number of elements of the next time level, which are
                                     adjacent to elements of the current time level.
   * \param[in]  adjElem           - The ID's of the adjacent elements.
   * \param[in]  secondPartTimeInt - Whether or not this is the second part of the
                                     time interval for the adjacent elements.
   */
  void ADER_DG_TimeInterpolatePredictorSol(CConfig             *config,
                                           const unsigned short iTime,
                                           const unsigned long  elemBeg,
                                           const unsigned long  elemEnd,
                                           const unsigned long  nAdjElem,
                                           const unsigned long  *adjElem,
                                           const bool           secondPartTimeInt);

  /*!
   * \brief Compute the artificial viscosity for shock capturing in DG. It is a virtual
            function, because this function is overruled for Navier-Stokes.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elemBeg - Begin index of the element range to be computed.
   * \param[in]  elemEnd - End index (not included) of the element range to be computed.
   */
  virtual void Shock_Capturing_DG(CConfig             *config,
                                  const unsigned long elemBeg,
                                  const unsigned long elemEnd);

  /*!
   * \brief Compute the volume contributions to the spatial residual. It is a virtual
            function, because this function is overruled for Navier-Stokes.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elemBeg - Begin index of the element range to be computed.
   * \param[in]  elemEnd - End index (not included) of the element range to be computed.
   */
  virtual void Volume_Residual(CConfig             *config,
                               const unsigned long elemBeg,
                               const unsigned long elemEnd);

  /*!
   * \brief Function, which computes the spatial residual for the DG discretization.
   * \param[in]  timeLevel           - Time level of the time accurate local time stepping,
                                       if relevant.
   * \param[in]  config              - Definition of the particular problem.
   * \param[in]  numerics            - Description of the numerical method.
   * \param[in]  haloInfoNeededForBC - If true,  treat boundaries for which halo data is needed.
                                       If false, treat boundaries for which only owned data is needed.
   */
  void Boundary_Conditions(const unsigned short timeLevel,
                           CConfig              *config,
                           CNumerics            **numerics,
                           const bool           haloInfoNeededForBC);

  /*!
   * \brief Compute the spatial residual for the given range of faces. It is a virtual
            function, because this function is overruled for Navier-Stokes.
   * \param[in] config      - Definition of the particular problem.
   * \param[in] indFaceBeg  - Starting index in the matching faces.
   * \param[in] indFaceEnd  - End index in the matching faces.
   * \param[in] numerics    - Description of the numerical method.
   */
  virtual void ResidualFaces(CConfig             *config,
                             const unsigned long indFaceBeg,
                             const unsigned long indFaceEnd,
                             CNumerics           *numerics);

  /*!
   * \brief Function, which accumulates the space time residual of the ADER-DG
            time integration scheme for the owned elements.
   * \param[in] config    - Definition of the particular problem.
   * \param[in] timeLevel - time level for which the residuals must be
                            accumulated.
   * \param[in] intPoint  - Index of the time integration point.
   */
  void AccumulateSpaceTimeResidualADEROwnedElem(CConfig             *config,
                                                const unsigned short timeLevel,
                                                const unsigned short intPoint);

  /*!
   * \brief Function, which accumulates the space time residual of the ADER-DG
            time integration scheme for the halo elements.
   * \param[in] config    - Definition of the particular problem.
   * \param[in] timeLevel - time level for which the residuals must be
                            accumulated.
   * \param[in] intPoint  - Index of the time integration point.
   */
  void AccumulateSpaceTimeResidualADERHaloElem(CConfig             *config,
                                               const unsigned short timeLevel,
                                               const unsigned short intPoint);

  /*!
   * \brief Compute primitive variables and their gradients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iStep - Current step in the time accurate local time
                        stepping algorithm, if appropriate.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void Preprocessing(CGeometry     *geometry,
                    CSolver        **solver_container,
                    CConfig        *config,
                    unsigned short iMesh,
                    unsigned short iRKStep,
                    unsigned short RunTime_EqSystem,
                    bool           Output) final;

  /*!
   * \brief
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Postprocessing(CGeometry      *geometry,
                      CSolver        **solver_container,
                      CConfig        *config,
                      unsigned short iMesh) final;

  /*!
   * \brief Impose via the residual the Euler wall boundary condition. It is a
            virtual function, because for Navier-Stokes it is overwritten.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   */
  virtual void BC_Euler_Wall(CConfig             *config,
                             const unsigned long surfElemBeg,
                             const unsigned long surfElemEnd,
                             CSurfaceElementFEM  *surfElem,
                             CNumerics           *conv_numerics);
  using CSolver::BC_Euler_Wall;

  /*!
   * \brief Impose the far-field boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   */
  virtual void BC_Far_Field(CConfig             *config,
                            const unsigned long surfElemBeg,
                            const unsigned long surfElemEnd,
                            CSurfaceElementFEM  *surfElem,
                            CNumerics           *conv_numerics);
  using CSolver::BC_Far_Field;

  /*!
   * \brief Impose the symmetry boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   */
  virtual void BC_Sym_Plane(CConfig             *config,
                            const unsigned long surfElemBeg,
                            const unsigned long surfElemEnd,
                            CSurfaceElementFEM  *surfElem,
                            CNumerics           *conv_numerics);
  using CSolver::BC_Sym_Plane;

  /*!
   * \brief Impose the supersonic outlet boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   */
  virtual void BC_Supersonic_Outlet(CConfig             *config,
                                    const unsigned long surfElemBeg,
                                    const unsigned long surfElemEnd,
                                    CSurfaceElementFEM  *surfElem,
                                    CNumerics           *conv_numerics);
  using CSolver::BC_Supersonic_Outlet;

  /*!
   * \brief Impose the subsonic inlet boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     val_marker    - Surface marker where the boundary condition is applied.
   */
  virtual void BC_Inlet(CConfig             *config,
                        const unsigned long surfElemBeg,
                        const unsigned long surfElemEnd,
                        CSurfaceElementFEM  *surfElem,
                        CNumerics           *conv_numerics,
                        unsigned short      val_marker);
  using CSolver::BC_Inlet;

  /*!
   * \brief Impose the outlet boundary condition.It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     val_marker    - Surface marker where the boundary condition is applied.
   */
  virtual void BC_Outlet(CConfig             *config,
                         const unsigned long surfElemBeg,
                         const unsigned long surfElemEnd,
                         CSurfaceElementFEM  *surfElem,
                         CNumerics           *conv_numerics,
                         unsigned short      val_marker);
  using CSolver::BC_Outlet;

  /*!
   * \brief Impose a constant heat-flux condition at the wall. It is a virtual
            function, such that it can be overwritten for Navier-Stokes.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     val_marker    - Surface marker where the boundary condition is applied.
   */
  virtual void BC_HeatFlux_Wall(CConfig            *config,
                                const unsigned long surfElemBeg,
                                const unsigned long surfElemEnd,
                                CSurfaceElementFEM *surfElem,
                                CNumerics           *conv_numerics,
                                unsigned short      val_marker);
  using CSolver::BC_HeatFlux_Wall;

  /*!
   * \brief Impose an isothermal condition at the wall. It is a virtual
            function, such that it can be overwritten for Navier-Stokes.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     val_marker    - Surface marker where the boundary condition is applied.
   */
  virtual void BC_Isothermal_Wall(CConfig            *config,
                                  const unsigned long surfElemBeg,
                                  const unsigned long surfElemEnd,
                                  CSurfaceElementFEM  *surfElem,
                                  CNumerics           *conv_numerics,
                                  unsigned short      val_marker);
  using CSolver::BC_Isothermal_Wall;

  /*!
   * \brief Impose the boundary condition using characteristic reconstruction. It is
   *        a virtual function, such that it can be overwritten for Navier-Stokes.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     val_marker    - Surface marker where the boundary condition is applied.
   */
  virtual void BC_Riemann(CConfig             *config,
                          const unsigned long surfElemBeg,
                          const unsigned long surfElemEnd,
                          CSurfaceElementFEM  *surfElem,
                          CNumerics           *conv_numerics,
                          unsigned short      val_marker);
  using CSolver::BC_Riemann;

  /*!
   * \brief Impose the user customized boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   */
  virtual void BC_Custom(CConfig             *config,
                         const unsigned long surfElemBeg,
                         const unsigned long surfElemEnd,
                         CSurfaceElementFEM  *surfElem,
                         CNumerics           *conv_numerics);
  using CSolver::BC_Custom;

  /*!
   * \brief Update the solution using a Runge-Kutta scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void ExplicitRK_Iteration(CGeometry *geometry,
                            CSolver **solver_container,
                            CConfig *config,
                            unsigned short iRKStep) final;

  /*!
   * \brief Update the solution using the classical fourth-order Runge-Kutta scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void ClassicalRK4_Iteration(CGeometry *geometry,
                              CSolver **solver_container,
                              CConfig *config,
                              unsigned short iRKStep) final;

  /*!
   * \brief Update the solution using the classical fourth-order Runge-Kutta scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetResidual_RMS_FEM(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Compute the global error measures (L2, Linf) for verification cases.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config   - Definition of the particular problem.
   */
  void ComputeVerificationError(CGeometry *geometry, CConfig *config) final;

  /*!
   * \brief Update the solution for the ADER-DG scheme for the given range
            of elements.
   * \param[in] elemBeg - Begin index of the element range to be computed.
   * \param[in] elemEnd - End index (not included) of the element range to be computed.
   */
  void ADER_DG_Iteration(const unsigned long elemBeg,
                         const unsigned long elemEnd);

  /*!
   * \brief Compute the pressure forces and all the adimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Pressure_Forces(const CGeometry* geometry, const CConfig* config) final;

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
                   bool val_update_geo) final;

  /*!
   * \brief The Euler and NS solvers support MPI+OpenMP.
   */
  inline bool GetHasHybridParallel() const final { return true; }

protected:

  /*!
   * \brief Routine that initiates the non-blocking communication between ranks
            for the givem time level.
   * \param[in] config    - Definition of the particular problem.
   * \param[in] timeLevel - The time level for which the communication must be
                            initiated.
   */
  void Initiate_MPI_Communication(CConfig *config,
                                  const unsigned short timeLevel);

  /*!
   * \brief Routine that initiates the reverse non-blocking communication
            between ranks.
   * \param[in] config    - Definition of the particular problem.
   * \param[in] timeLevel - The time level for which the reverse communication
                            must be initiated.
   */
  void Initiate_MPI_ReverseCommunication(CConfig *config,
                                         const unsigned short timeLevel);

  /*!
   * \brief Routine that completes the non-blocking communication between ranks.
   * \param[in] config              - Definition of the particular problem.
   * \param[in] timeLevel           - The time level for which the communication
                                      may be completed.
   * \param[in] commMustBeCompleted - Whether or not the communication must be completed.
   * \return  Whether or not the communication has been completed.
   */
  bool Complete_MPI_Communication(CConfig *config,
                                  const unsigned short timeLevel,
                                  const bool commMustBeCompleted);

  /*!
   * \brief Routine that completes the reverse non-blocking communication
            between ranks.
   * \param[in] config              - Definition of the particular problem.
   * \param[in] timeLevel           - The time level for which the communication
                                      may be completed.
   * \param[in] commMustBeCompleted - Whether or not the communication must be completed.
   * \return  Whether or not the communication has been completed.
   */
  bool Complete_MPI_ReverseCommunication(CConfig *config,
                                         const unsigned short timeLevel,
                                         const bool commMustBeCompleted);

  /*!
   * \brief Function, which computes the inviscid fluxes in face points.
   * \param[in]  config         - Definition of the particular problem.
   * \param[in]  solLeft        - Left solution in the points.
   * \param[in]  solRight       - Right solution in the points.
   * \param[in]  JacobiansFace  - The jacobians of the face.
   * \param[in]  normalsFace    - The normals of the face.
   * \param[in]  gridVelocities - The grid velocities of the face.
   * \param[in]  numerics       - Object, which contains the Riemann solver.
   * \param[out] fluxes         - The fluxes to be computed.
   */
  void ComputeInviscidFluxesFace(CConfig                   *config,
                                 ColMajorMatrix<su2double> &solLeft,
                                 ColMajorMatrix<su2double> &solRight,
                                 su2activevector           &JacobiansFace,
                                 ColMajorMatrix<su2double> &normalsFace,
                                 ColMajorMatrix<su2double> &gridVelocities,
                                 CNumerics                 *numerics,
                                 ColMajorMatrix<su2double> &fluxes);

  /*!
   * \brief Function, which computes the left state of a boundary face.
   * \param[in]  config     - Definition of the particular problem.
   * \param[in]  nFaceSimul - Number of faces that are treated simultaneously
                              to improve performance.
   * \param[in]  NPad       - Value of the padding parameter to obtain optimal
                              performance in the gemm computations.
   * \param[in]  surfElem   - Surface boundary elements for which the left state must be computed.
   * \param[out] solFace    - Temporary storage for the solution in the DOFs.
   * \param[out] solIntL    - Left states in the integration points of the face.
   */
  void LeftStatesIntegrationPointsBoundaryFace(CConfig                  *config,
                                               const unsigned short     nFaceSimul,
                                               const unsigned short     NPad,
                                               const CSurfaceElementFEM *surfElem,
                                               su2double                *solFace,
                                               su2double                *solIntL);

  /*!
   * \brief Function, which computes the boundary states for the given items
            of the boundary face by applying the inviscid wall boundary conditions.
   * \param[in]  surfElem   - Surface boundary elements for which the left state must
                              be computed.
   * \param[in]  solIntL    - Left states in the integration points of the face.
   * \param[out] solIntR    - Right states in the integration points of the face.
   */
  void BoundaryStates_Euler_Wall(const CSurfaceElementFEM        *surfElem,
                                 const ColMajorMatrix<su2double> &solIntL,
                                 ColMajorMatrix<su2double>       &solIntR);

  /*!
   * \brief Function, which computes the boundary states in the integration points
            of the boundary face by applying the inlet boundary conditions.
   * \param[in]  config     - Definition of the particular problem.
   * \param[in]  nFaceSimul - Number of fused faces, i.e. the number of faces
                              that are treated simultaneously to improve performance.
   * \param[in]  NPad       - Value of the padding parameter to obtain optimal
                              performance in the gemm computations.
   * \param[in]  surfElem   - Surface boundary element for which the left state must be computed.
   * \param[in]  val_marker - Surface marker where the boundary condition is applied.
   * \param[in]  solIntL    - Left states in the integration points of the face.
   * \param[out] solIntR    - Right states in the integration points of the face.
   */
  void BoundaryStates_Inlet(CConfig                  *config,
                            const unsigned short     nFaceSimul,
                            const unsigned short     NPad,
                            const CSurfaceElementFEM *surfElem,
                            unsigned short           val_marker,
                            const su2double          *solIntL,
                            su2double                *solIntR);

  /*!
   * \brief Function, which computes the boundary states in the integration points
            of the boundary face by applying the outlet boundary conditions.
   * \param[in]  config     - Definition of the particular problem.
   * \param[in]  nFaceSimul - Number of fused faces, i.e. the number of faces
                              that are treated simultaneously to improve performance.
   * \param[in]  NPad       - Value of the padding parameter to obtain optimal
                              performance in the gemm computations.
   * \param[in]  surfElem   - Surface boundary element for which the left state must be computed.
   * \param[in]  val_marker - Surface marker where the boundary condition is applied.
   * \param[in]  solIntL    - Left states in the integration points of the face.
   * \param[out] solIntR    - Right states in the integration points of the face.
   */
  void BoundaryStates_Outlet(CConfig                  *config,
                             const unsigned short     nFaceSimul,
                             const unsigned short     NPad,
                             const CSurfaceElementFEM *surfElem,
                             unsigned short           val_marker,
                             const su2double          *solIntL,
                             su2double                *solIntR);

  /*!
   * \brief Function, which computes the boundary states in the integration points
            of the boundary face by applying the Riemann boundary conditions.
   * \param[in]  config     - Definition of the particular problem.
   * \param[in]  nFaceSimul - Number of fused faces, i.e. the number of faces
                              that are treated simultaneously to improve performance.
   * \param[in]  NPad       - Value of the padding parameter to obtain optimal
                              performance in the gemm computations.
   * \param[in]  surfElem   - Surface boundary element for which the left state must be computed.
   * \param[in]  val_marker - Surface marker where the boundary condition is applied.
   * \param[in]  solIntL    - Left states in the integration points of the face.
   * \param[out] solIntR    - Right states in the integration points of the face.
   */
  void BoundaryStates_Riemann(CConfig                  *config,
                              const unsigned short     nFaceSimul,
                              const unsigned short     NPad,
                              const CSurfaceElementFEM *surfElem,
                              unsigned short           val_marker,
                              const su2double          *solIntL,
                              su2double                *solIntR);

  /*!
   * \brief Function, which converts the entropy variables to the
   *        primitive variables for the given entities.
   * \param[in,out] sol  - Matrix, which contains the entropy variables
   *                       on entry and the primitive variables on exit.
   */
  void EntropyToPrimitiveVariables(ColMajorMatrix<su2double> &sol);

  /*!
   * \brief Function, which computes the volume source terms for the given element.
   * \param[in]     config - Definition of the particular problem.
   * \param[in]     elem   - Volume element for which the source terms must be computed.
   * \param[in,out] sol    - Matrix, which contains the primitive variables in the
   *                         integration points on entry and the source terms on exit,
   *                         if these must be computed.
   * \return - Whether or not source terms are present.
   */
  bool VolumeSourceTerms(CConfig                  *config,
                         CVolumeElementFEM_DG     *elem,
                         ColMajorMatrix<su2double> &sol);

  /*!
   * \brief Function, which computes the transformation matrix between conservative
   *        and entropy variables in the integration points of the element.
   * \param[in,out] elem - Volume element for which the transformation matrix
   *                       must be computed.
   * \param[in] sol     - Matrix, which contains the primitive variables in the
   *                      integration points.
   */
  void VolumeTransformationMatrix(CVolumeElementFEM_DG     *elem,
                                  ColMajorMatrix<su2double> &sol);
private:

  /*!
   * \brief Virtual function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using an
            aliased discretization in 2D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  virtual void ADER_DG_AliasedPredictorResidual_2D(CConfig              *config,
                                                   CVolumeElementFEM_DG *elem,
                                                   const su2double      *sol,
                                                   const unsigned short nSimul,
                                                   const unsigned short NPad,
                                                   su2double            *res,
                                                   su2double            *work);

  /*!
   * \brief Virtual function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using an
            aliased discretization in 3D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  virtual void ADER_DG_AliasedPredictorResidual_3D(CConfig              *config,
                                                   CVolumeElementFEM_DG *elem,
                                                   const su2double      *sol,
                                                   const unsigned short nSimul,
                                                   const unsigned short NPad,
                                                   su2double            *res,
                                                   su2double            *work);

  /*!
   * \brief Virtual function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using a
            non-aliased discretization in 2D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  virtual void ADER_DG_NonAliasedPredictorResidual_2D(CConfig              *config,
                                                      CVolumeElementFEM_DG *elem,
                                                      const su2double      *sol,
                                                      const unsigned short nSimul,
                                                      const unsigned short NPad,
                                                      su2double            *res,
                                                      su2double            *work);

  /*!
   * \brief Virtual function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using a
            non-aliased discretization in 3D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  virtual void ADER_DG_NonAliasedPredictorResidual_3D(CConfig              *config,
                                                      CVolumeElementFEM_DG *elem,
                                                      const su2double      *sol,
                                                      const unsigned short nSimul,
                                                      const unsigned short NPad,
                                                      su2double            *res,
                                                      su2double            *work);

  /*!
   * \brief Function, which converts the conservative variables to the
   *        entropy variables for the given entities.
   * \param[in,out] sol  - Matrix, which contains the conservative variables
   *                       on entry and the entropy variables on exit.
   */
  void ConservativeToEntropyVariables(ColMajorMatrix<su2double> &sol);

  /*!
   * \brief Function, which sets up the list of tasks to be carried out in the
            computationally expensive part of the solver.
   * \param[in] config - Definition of the particular problem.
   */
  void SetUpTaskList(CConfig *config);

  /*!
   * \brief Function, which creates the final residual by summing up
            the contributions for the DOFs of the elements considered.
   * \param[in] timeLevel     - Time level of the elements for which the
                                final residual must be created.
   * \param[in] ownedElements - Whether owned or halo elements must be treated.
   */
  void CreateFinalResidual(const unsigned short timeLevel,
                           const bool ownedElements);

  /*!
   * \brief Function, which multiplies the residual by the inverse
            of the (lumped) mass matrix.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  useADER - Whether or not the ADER residual must be multiplied.
   * \param[in]  elemBeg - Begin index of the element range to be computed.
   * \param[in]  elemEnd - End index (not included) of the element range to be computed.
   */
  void MultiplyResidualByInverseMassMatrix(CConfig             *config,
                                           const bool          useADER,
                                           const unsigned long elemBeg,
                                           const unsigned long elemEnd);

  /*!
   * \brief Function, which computes the residual contribution from a boundary
            face in an inviscid computation when the boundary conditions have
            already been applied.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     surfElem      - Surface boundary element for which the
                                    contribution to the residual must be computed.
   * \param[in]     solInt0       - Solution in the integration points of side 0.
   * \param[in]     solInt1       - Solution in the integration points of side 1.
   */
  void ResidualInviscidBoundaryFace(CConfig                   *config,
                                    CNumerics                 *conv_numerics,
                                    CSurfaceElementFEM        *surfElem,
                                    ColMajorMatrix<su2double> &solInt0,
                                    ColMajorMatrix<su2double> &solInt1);
};
