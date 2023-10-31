/*!
 * \file CFEM_DG_EulerSolver.hpp
 * \brief Headers of the CFEM_DG_EulerSolver class
 * \author E. van der Weide, T. Economon, J. Alonso
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

/*!
 * \class CFEM_DG_EulerSolver
 * \brief Main class for defining the Euler Discontinuous Galerkin finite element flow solver.
 * \ingroup Euler_Equations
 * \author E. van der Weide, T. Economon, J. Alonso
 * \version 8.0.0 "Harrier"
 */
class CFEM_DG_EulerSolver : public CSolver {
protected:

  su2double Gamma;           /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One; /*!< \brief Fluids's Gamma - 1.0  . */

  CFluidModel  *FluidModel; /*!< \brief fluid model used in the solver */

  su2double
  Mach_Inf,         /*!< \brief Mach number at infinity. */
  Density_Inf,      /*!< \brief Density at infinity. */
  Energy_Inf,     /*!< \brief Energy at infinity. */
  Temperature_Inf,  /*!< \brief Energy at infinity. */
  Pressure_Inf,     /*!< \brief Pressure at infinity. */
  *Velocity_Inf;    /*!< \brief Flow velocity vector at infinity. */

  vector<su2double> ConsVarFreeStream; /*!< \brief Vector, which contains the free stream
                                        conservative variables. */
  su2double
  *CL_Inv,        /*!< \brief Lift coefficient (inviscid contribution) for each boundary. */
  *CD_Inv,        /*!< \brief Drag coefficient (inviscid contribution) for each boundary. */
  *CSF_Inv,       /*!< \brief Sideforce coefficient (inviscid contribution) for each boundary. */
  *CFx_Inv,       /*!< \brief x Force coefficient (inviscid contribution) for each boundary. */
  *CFy_Inv,       /*!< \brief y Force coefficient (inviscid contribution) for each boundary. */
  *CFz_Inv,       /*!< \brief z Force coefficient (inviscid contribution) for each boundary. */
  *CMx_Inv,       /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CMy_Inv,       /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CMz_Inv,       /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CEff_Inv;      /*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each boundary. */

  su2double
  *Surface_CL_Inv,    /*!< \brief Lift coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CD_Inv,    /*!< \brief Drag coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CSF_Inv,   /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFx_Inv,   /*!< \brief x Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFy_Inv,   /*!< \brief y Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFz_Inv,   /*!< \brief z Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMx_Inv,   /*!< \brief x Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMy_Inv,   /*!< \brief y Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMz_Inv,   /*!< \brief z Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CEff_Inv;  /*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each monitoring surface. */

  su2double
  AllBound_CL_Inv,    /*!< \brief Total lift coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CD_Inv,    /*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CSF_Inv,   /*!< \brief Total sideforce coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFx_Inv,   /*!< \brief Total x force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Inv,   /*!< \brief Total y force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Inv,   /*!< \brief Total z force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMx_Inv,   /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Inv,   /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Inv,   /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Inv;  /*!< \brief Total efficiency (Cl/Cd) (inviscid contribution) for all the boundaries. */

  su2double
  AeroCoeffForceRef,   /*!< \brief Reference force for coefficients */
  Total_CL,      /*!< \brief Total lift coefficient for all the boundaries. */
  Total_CD,        /*!< \brief Total drag coefficient for all the boundaries. */
  Total_CSF,       /*!< \brief Total sideforce coefficient for all the boundaries. */
  Total_CFx,     /*!< \brief Total x force coefficient for all the boundaries. */
  Total_CFy,     /*!< \brief Total y force coefficient for all the boundaries. */
  Total_CFz,     /*!< \brief Total z force coefficient for all the boundaries. */
  Total_CMx,     /*!< \brief Total x moment coefficient for all the boundaries. */
  Total_CMy,     /*!< \brief Total y moment coefficient for all the boundaries. */
  Total_CMz,     /*!< \brief Total z moment coefficient for all the boundaries. */
  Total_CEff;      /*!< \brief Total efficiency coefficient for all the boundaries. */

  su2double
  *Surface_CL,         /*!< \brief Lift coefficient for each monitoring surface. */
  *Surface_CD,         /*!< \brief Drag coefficient for each monitoring surface. */
  *Surface_CSF,        /*!< \brief Side-force coefficient for each monitoring surface. */
  *Surface_CFx,        /*!< \brief x Force coefficient for each monitoring surface. */
  *Surface_CFy,        /*!< \brief y Force coefficient for each monitoring surface. */
  *Surface_CFz,        /*!< \brief z Force coefficient for each monitoring surface. */
  *Surface_CMx,        /*!< \brief x Moment coefficient for each monitoring surface. */
  *Surface_CMy,        /*!< \brief y Moment coefficient for each monitoring surface. */
  *Surface_CMz,        /*!< \brief z Moment coefficient for each monitoring surface. */
  *Surface_CEff;       /*!< \brief Efficiency (Cl/Cd) for each monitoring surface. */

  unsigned long nDOFsLocTot;    /*!< \brief Total number of local DOFs, including halos. */
  unsigned long nDOFsLocOwned;  /*!< \brief Number of owned local DOFs. */
  unsigned long nDOFsGlobal;    /*!< \brief Number of global DOFs. */

  unsigned long nVolElemTot;    /*!< \brief Total number of local volume elements, including halos. */
  unsigned long nVolElemOwned;  /*!< \brief Number of owned local volume elements. */
  CVolumeElementFEM *volElem;   /*!< \brief Array of the local volume elements, including halos. */

  const unsigned long *nVolElemOwnedPerTimeLevel;    /*!< \brief Number of owned local volume elements
                                                                 per time level. Cumulative storage. */
  const unsigned long *nVolElemInternalPerTimeLevel; /*!< \brief Number of internal local volume elements per
                                                                 time level. Internal means that the solution
                                                                 data does not need to be communicated. */
  const unsigned long *nVolElemHaloPerTimeLevel;     /*!< \brief Number of halo volume elements
                                                                 per time level. Cumulative storage. */

  vector<vector<unsigned long> > ownedElemAdjLowTimeLevel; /*!< \brief List of owned elements per time level that are
                                                                       adjacent to elements of the lower time level. */
  vector<vector<unsigned long> > haloElemAdjLowTimeLevel;  /*!< \brief List of halo elements per time level that are
                                                                       adjacent to elements of the lower time level. */

  unsigned long nMeshPoints;    /*!< \brief Number of mesh points in the local part of the grid. */
  CPointFEM *meshPoints;        /*!< \brief Array of the points of the FEM mesh. */

  const unsigned long *nMatchingInternalFacesWithHaloElem;  /*!< \brief Number of local matching internal faces per time level
                                                                        between an owned and a halo element. Cumulative storage. */
  const unsigned long *nMatchingInternalFacesLocalElem;     /*!< \brief Number of local matching internal faces per time level
                                                                        between local elements. Cumulative storage. */

  CInternalFaceElementFEM *matchingInternalFaces;    /*!< \brief Array of the local matching internal faces. */
  CBoundaryFEM *boundaries;                          /*!< \brief Array of the boundaries of the FEM mesh. */

  unsigned short nStandardBoundaryFacesSol; /*!< \brief Number of standard boundary faces used for solution of the DG solver. */
  unsigned short nStandardElementsSol;      /*!< \brief Number of standard volume elements used for solution of the DG solver. */
  unsigned short nStandardMatchingFacesSol; /*!< \brief Number of standard matching internal faces used for solution of the DG solver. */

  const CFEMStandardBoundaryFace *standardBoundaryFacesSol; /*!< \brief Array that contains the standard boundary
                                                                        faces used for the solution of the DG solver. */
  const CFEMStandardElement *standardElementsSol; /*!< \brief Array that contains the standard volume elements
                                                              used for the solution of the DG solver. */
  const CFEMStandardInternalFace *standardMatchingFacesSol;  /*!< \brief Array that contains the standard matching
                                                                         internal faces used for the solution of
                                                                         the DG solver. */

  const su2double *timeCoefADER_DG;                        /*!< \brief The time coefficients in the iteration matrix of
                                                                       the ADER-DG predictor step. */
  const su2double *timeInterpolDOFToIntegrationADER_DG;    /*!< \brief The interpolation matrix between the time DOFs and
                                                                       the time integration points for ADER-DG. */
  const su2double *timeInterpolAdjDOFToIntegrationADER_DG; /*!< \brief The interpolation matrix between the time DOFs of adjacent
                                                                       elements of a higher time level and the time integration
                                                                       points for ADER-DG. */

  unsigned int sizeWorkArray;     /*!< \brief The size of the work array needed. */

  vector<su2double> TolSolADER;   /*!< \brief Vector, which stores the tolerances for the conserved
                                              variables in the ADER predictor step. */

  vector<su2double> VecSolDOFs;    /*!< \brief Vector, which stores the solution variables in the owned DOFs. */
  vector<su2double> VecSolDOFsNew; /*!< \brief Vector, which stores the new solution variables in the owned DOFs (needed for classical RK4 scheme). */
  vector<su2double> VecDeltaTime;  /*!< \brief Vector, which stores the time steps of the owned volume elements. */

  vector<su2double> VecSolDOFsPredictorADER; /*!< \brief Vector, which stores the ADER predictor solution in the owned
                                                         DOFs. These are both space and time DOFs. */

  vector<vector<su2double> > VecWorkSolDOFs; /*!< \brief Working double vector to store the conserved variables for
                                                         the DOFs for the different time levels. */

  vector<su2double> VecResDOFs;         /*!< \brief Vector, which stores the residuals in the owned DOFs. */
  vector<su2double> VecResFaces;        /*!< \brief Vector, which stores the residuals of the DOFs that
                                                    come from the faces, both boundary and internal. */
  vector<su2double> VecTotResDOFsADER;  /*!< \brief Vector, which stores the accumulated residuals of the
                                                    owned DOFs for the ADER corrector step. */


  vector<unsigned long> nEntriesResFaces; /*!< \brief Number of entries for the DOFs in the
                                                      residual of the faces. Cumulative storage. */
  vector<unsigned long> entriesResFaces;  /*!< \brief The corresponding entries in the residual of the faces. */

  vector<unsigned long> nEntriesResAdjFaces; /*!< \brief Number of entries for the DOFs in the residual of the faces,
                                                         where the face is adjacent to an element of lower time
                                                         level. Cumulative storage. */
  vector<unsigned long> entriesResAdjFaces;  /*!< \brief The corresponding entries in the residual of the faces. */

  vector<vector<unsigned long> > startLocResFacesMarkers; /*!< \brief The starting location in the residual of the
                                                                      faces for the time levels of the boundary
                                                                      markers. */

  vector<unsigned long> startLocResInternalFacesLocalElem; /*!< \brief The starting location in the residual of the
                                                                       faces for the time levels of internal faces
                                                                       between locally owned elements. */
  vector<unsigned long> startLocResInternalFacesWithHaloElem; /*!< \brief The starting location in the residual of the
                                                                          faces for the time levels of internal faces
                                                                          between an owned and a halo element. */

  bool symmetrizingTermsPresent;    /*!< \brief Whether or not symmetrizing terms are present in the
                                                discretization. */

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

  CBlasStructure *blasFunctions; /*!< \brief  Pointer to the object to carry out the BLAS functionalities. */

private:

#ifdef HAVE_MPI
  vector<vector<SU2_MPI::Request> > commRequests;  /*!< \brief Communication requests in the communication of the solution for all
                                                               time levels. These are both sending and receiving requests. */

  vector<vector<vector<unsigned long> > > elementsRecvMPIComm;  /*!< \brief Triple vector, which contains the halo elements
                                                                            for MPI communication for all time levels. */
  vector<vector<vector<unsigned long> > > elementsSendMPIComm;  /*!< \brief Triple vector, which contains the donor elements
                                                                            for MPI communication for all time levels. */

  vector<vector<int> > ranksRecvMPI; /*!< \brief Double vector, which contains the ranks from which the halo elements
                                                 are received for all time levels. */
  vector<vector<int> > ranksSendMPI; /*!< \brief Double vector, which contains the ranks to which the donor elements
                                                 are sent for all time levels. */

  vector<vector<vector<su2double> > > commRecvBuf;  /*!< \brief Receive buffers used to receive the solution data
                                                                in the communication pattern for all time levels. */
  vector<vector<vector<su2double> > > commSendBuf;  /*!< \brief Send buffers used to send the solution data
                                                                in the communication pattern for all time levels. */
#endif

  vector<vector<unsigned long> > elementsRecvSelfComm;  /*!< \brief Double vector, which contains the halo elements
                                                                    for self communication for all time levels. */
  vector<vector<unsigned long> > elementsSendSelfComm;  /*!< \brief Double vector, which contains the donor elements
                                                                    for self communication for all time levels. */

  vector<su2double> rotationMatricesPeriodicity;    /*!< \brief Vector, which contains the rotation matrices
                                                                for the rotational periodic transformations. */
  vector<vector<vector<unsigned long> > > halosRotationalPeriodicity; /*!< \brief Triple vector, which contains the indices
                                                                                  of halo elements for which a periodic
                                                                                  transformation must be applied for all
                                                                                  time levels. */

  vector<CTaskDefinition> tasksList; /*!< \brief List of tasks to be carried out in the computationally
                                                 intensive part of the solver. */

  CVariable* GetBaseClassPointerToNodes() final {return nullptr;}

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
   * \brief Get a pointer to the vector of the solution degrees of freedom.
   * \return Pointer to the vector of the solution degrees of freedom.
   */
  inline su2double* GetVecSolDOFs(void) final { return VecSolDOFs.data(); }

  /*!
   * \brief Get the global number of solution degrees of freedom for the calculation.
   * \return Global number of solution degrees of freedom
   */
  inline unsigned long GetnDOFsGlobal(void) const final { return nDOFsGlobal; }

  /*!
   * \brief Compute the pressure at the infinity.
   * \return Value of the pressure at the infinity.
   */
  inline CFluidModel* GetFluidModel(void) const final { return FluidModel;}

  /*!
   * \brief Compute the density at the infinity.
   * \return Value of the density at the infinity.
   */
  inline su2double GetDensity_Inf(void) const final { return Density_Inf; }

  /*!
   * \brief Compute 2-norm of the velocity at the infinity.
   * \return Value of the 2-norm of the velocity at the infinity.
   */
  inline su2double GetModVelocity_Inf(void) const final {
    su2double Vel2 = 0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    return sqrt(Vel2);
  }


  /*!
   * \brief Compute the density multiply by energy at the infinity.
   * \return Value of the density multiply by  energy at the infinity.
   */
  inline su2double GetDensity_Energy_Inf(void) const final { return Density_Inf*Energy_Inf; }

  /*!
   * \brief Compute the pressure at the infinity.
   * \return Value of the pressure at the infinity.
   */
  inline su2double GetPressure_Inf(void) const final { return Pressure_Inf; }

  /*!
   * \brief Get the velocity at the infinity.
   * \param[in] val_dim - Index of the velocity vector.
   * \return Value of the velocity at the infinity.
   */
  inline su2double GetVelocity_Inf(unsigned short val_dim) const final { return Velocity_Inf[val_dim]; }

  /*!
   * \brief Get the velocity at the infinity.
   * \return Value of the velocity at the infinity.
   */
  inline su2double *GetVelocity_Inf(void) final { return Velocity_Inf; }

  /*!
   * \brief Set the freestream pressure.
   * \param[in] Value of freestream pressure.
   */
  inline void SetPressure_Inf(su2double p_inf) final { Pressure_Inf = p_inf; }

  /*!
   * \brief Set the freestream temperature.
   * \param[in] Value of freestream temperature.
   */
  inline void SetTemperature_Inf(su2double t_inf) final { Temperature_Inf = t_inf; }

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
   * \param[out] workArray - Work array.
   */
  void ADER_DG_PredictorStep(CConfig             *config,
                             const unsigned long elemBeg,
                             const unsigned long elemEnd,
                             su2double           *workArray);

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
   * \param[out] solTimeLevel      - Array in which the interpolated solution for the
                                     time level considered must be stored.
   */
  void ADER_DG_TimeInterpolatePredictorSol(CConfig             *config,
                                           const unsigned short iTime,
                                           const unsigned long  elemBeg,
                                           const unsigned long  elemEnd,
                                           const unsigned long  nAdjElem,
                                           const unsigned long  *adjElem,
                                           const bool           secondPartTimeInt,
                                           su2double            *solTimeLevel);

  /*!
   * \brief Compute the artificial viscosity for shock capturing in DG. It is a virtual
            function, because this function is overruled for Navier-Stokes.
   * \param[in]  config    - Definition of the particular problem.
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   * \param[out] workArray - Work array.
   */
  virtual void Shock_Capturing_DG(CConfig             *config,
                                  const unsigned long elemBeg,
                                  const unsigned long elemEnd,
                                  su2double           *workArray);

  /*!
   * \brief Compute the volume contributions to the spatial residual. It is a virtual
            function, because this function is overruled for Navier-Stokes.
   * \param[in]  config    - Definition of the particular problem.
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   * \param[out] workArray - Work array.
   */
  virtual void Volume_Residual(CConfig             *config,
                               const unsigned long elemBeg,
                               const unsigned long elemEnd,
                               su2double           *workArray);

  /*!
   * \brief Function, which computes the spatial residual for the DG discretization.
   * \param[in]  timeLevel           - Time level of the time accurate local time stepping,
                                       if relevant.
   * \param[in]  config              - Definition of the particular problem.
   * \param[in]  numerics            - Description of the numerical method.
   * \param[in]  haloInfoNeededForBC - If true,  treat boundaries for which halo data is needed.
                                       If false, treat boundaries for which only owned data is needed.
   * \param[out] workArray           - Work array.
   */
  void Boundary_Conditions(const unsigned short timeLevel,
                           CConfig              *config,
                           CNumerics            **numerics,
                           const bool           haloInfoNeededForBC,
                           su2double            *workArray);

  /*!
   * \brief Compute the spatial residual for the given range of faces. It is a virtual
            function, because this function is overruled for Navier-Stokes.
   * \param[in]     config      - Definition of the particular problem.
   * \param[in]     indFaceBeg  - Starting index in the matching faces.
   * \param[in]     indFaceEnd  - End index in the matching faces.
   * \param[in,out] indResFaces - Index where to store the residuals in
                                  the vector of face residuals.
   * \param[in]     numerics    - Description of the numerical method.
   * \param[out]    workArray   - Work array.
   */
  virtual void ResidualFaces(CConfig             *config,
                             const unsigned long indFaceBeg,
                             const unsigned long indFaceEnd,
                             unsigned long       &indResFaces,
                             CNumerics           *numerics,
                             su2double           *workArray);

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
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Euler_Wall(CConfig                  *config,
                             const unsigned long      surfElemBeg,
                             const unsigned long      surfElemEnd,
                             const CSurfaceElementFEM *surfElem,
                             su2double                *resFaces,
                             CNumerics                *conv_numerics,
                             su2double                *workArray);
  using CSolver::BC_Euler_Wall;

  /*!
   * \brief Impose the far-field boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Far_Field(CConfig                  *config,
                            const unsigned long      surfElemBeg,
                            const unsigned long      surfElemEnd,
                            const CSurfaceElementFEM *surfElem,
                            su2double                *resFaces,
                            CNumerics                *conv_numerics,
                            su2double                *workArray);
  using CSolver::BC_Far_Field;

  /*!
   * \brief Impose the symmetry boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Sym_Plane(CConfig                  *config,
                            const unsigned long      surfElemBeg,
                            const unsigned long      surfElemEnd,
                            const CSurfaceElementFEM *surfElem,
                            su2double                *resFaces,
                            CNumerics                *conv_numerics,
                            su2double                *workArray);
  using CSolver::BC_Sym_Plane;

  /*!
   * \brief Impose the supersonic outlet boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Supersonic_Outlet(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics,
                                    su2double                *workArray);
  using CSolver::BC_Supersonic_Outlet;

  /*!
   * \brief Impose the subsonic inlet boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Inlet(CConfig                  *config,
                        const unsigned long      surfElemBeg,
                        const unsigned long      surfElemEnd,
                        const CSurfaceElementFEM *surfElem,
                        su2double                *resFaces,
                        CNumerics                *conv_numerics,
                        unsigned short           val_marker,
                        su2double                *workArray);
  using CSolver::BC_Inlet;

  /*!
   * \brief Impose the outlet boundary condition.It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Outlet(CConfig                  *config,
                         const unsigned long      surfElemBeg,
                         const unsigned long      surfElemEnd,
                         const CSurfaceElementFEM *surfElem,
                         su2double                *resFaces,
                         CNumerics                *conv_numerics,
                         unsigned short           val_marker,
                         su2double                *workArray);
  using CSolver::BC_Outlet;

  /*!
   * \brief Impose a constant heat-flux condition at the wall. It is a virtual
            function, such that it can be overwritten for Navier-Stokes.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  inline virtual void BC_HeatFlux_Wall(CConfig *config,
                                      const unsigned long      surfElemBeg,
                                      const unsigned long      surfElemEnd,
                                      const CSurfaceElementFEM *surfElem,
                                      su2double                *resFaces,
                                      CNumerics                *conv_numerics,
                                      unsigned short           val_marker,
                                      su2double                *workArray) {}
  using CSolver::BC_HeatFlux_Wall;

  /*!
   * \brief Impose an isothermal condition at the wall. It is a virtual
            function, such that it can be overwritten for Navier-Stokes.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  inline virtual void BC_Isothermal_Wall(CConfig *config,
                                        const unsigned long      surfElemBeg,
                                        const unsigned long      surfElemEnd,
                                        const CSurfaceElementFEM *surfElem,
                                        su2double                *resFaces,
                                        CNumerics                *conv_numerics,
                                        unsigned short           val_marker,
                                        su2double                *workArray) {}
  using CSolver::BC_Isothermal_Wall;

  /*!
   * \brief Impose the boundary condition using characteristic reconstruction. It is
   *        a virtual function, such that it can be overwritten for Navier-Stokes.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Riemann(CConfig                  *config,
                          const unsigned long      surfElemBeg,
                          const unsigned long      surfElemEnd,
                          const CSurfaceElementFEM *surfElem,
                          su2double                *resFaces,
                          CNumerics                *conv_numerics,
                          unsigned short           val_marker,
                          su2double                *workArray);
  using CSolver::BC_Riemann;

  /*!
   * \brief Impose the user customized boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Custom(CConfig                  *config,
                         const unsigned long      surfElemBeg,
                         const unsigned long      surfElemEnd,
                         const CSurfaceElementFEM *surfElem,
                         su2double                *resFaces,
                         CNumerics                *conv_numerics,
                         su2double                *workArray);
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
   * \brief Provide the non dimensional lift coefficient (inviscid contribution).
   * \param val_marker Surface where the coefficient is going to be computed.
   * \return Value of the lift coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCL_Inv(unsigned short val_marker) const final { return CL_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CL(unsigned short val_marker) const final { return Surface_CL[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CD(unsigned short val_marker) const final { return Surface_CD[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CSF(unsigned short val_marker) const final { return Surface_CSF[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CEff(unsigned short val_marker) const final { return Surface_CEff[val_marker]; }

  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFx(unsigned short val_marker) const final { return Surface_CFx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFy(unsigned short val_marker) const final { return Surface_CFy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFz(unsigned short val_marker) const final { return Surface_CFz[val_marker]; }

  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMx(unsigned short val_marker) const final { return Surface_CMx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMy(unsigned short val_marker) const final { return Surface_CMy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMz(unsigned short val_marker) const final { return Surface_CMz[val_marker]; }

  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CL_Inv(unsigned short val_marker) const final { return Surface_CL_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CD_Inv(unsigned short val_marker) const final { return Surface_CD_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CSF_Inv(unsigned short val_marker) const final { return Surface_CSF_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CEff_Inv(unsigned short val_marker) const final { return Surface_CEff_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFx_Inv(unsigned short val_marker) const final { return Surface_CFx_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFy_Inv(unsigned short val_marker) const final { return Surface_CFy_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFz_Inv(unsigned short val_marker) const final { return Surface_CFz_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMx_Inv(unsigned short val_marker) const final { return Surface_CMx_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMy_Inv(unsigned short val_marker) const final { return Surface_CMy_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMz_Inv(unsigned short val_marker) const final { return Surface_CMz_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the drag coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCD_Inv(unsigned short val_marker) const final { return CD_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional sideforce coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the sideforce coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCSF_Inv(unsigned short val_marker) const final { return CSF_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional efficiency coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the efficiency coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCEff_Inv(unsigned short val_marker) const final { return CEff_Inv[val_marker]; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CSF() const final { return Total_CSF; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CEff() const final { return Total_CEff; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional lift coefficient.
   * \param[in] val_Total_CL - Value of the total lift coefficient.
   */
  inline void SetTotal_CL(su2double val_Total_CL) final { Total_CL = val_Total_CL; }

  /*!
   * \brief Get the reference force used to compute CL, CD, etc.
   */
  inline su2double GetAeroCoeffsReferenceForce() const final { return AeroCoeffForceRef; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional lift coefficient.
   * \return Value of the lift coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CL() const final { return Total_CL; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CD() const final { return Total_CD; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional x moment coefficient.
   * \return Value of the moment x coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMx() const final { return Total_CMx; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional y moment coefficient.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMy() const final { return Total_CMy; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional z moment coefficient.
   * \return Value of the moment z coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMz() const final { return Total_CMz; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional x force coefficient.
   * \return Value of the force x coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CFx() const final { return Total_CFx; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional y force coefficient.
   * \return Value of the force y coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CFy() const final { return Total_CFy; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional z force coefficient.
   * \return Value of the force z coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CFz() const final { return Total_CFz; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_CD(su2double val_Total_CD) final { Total_CD = val_Total_CD; }

  /*!
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CL_Inv() const final { return AllBound_CL_Inv; }

  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CD_Inv() const final { return AllBound_CD_Inv; }

  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CSF_Inv() const final { return AllBound_CSF_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CEff_Inv() const final { return AllBound_CEff_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMx_Inv() const final { return AllBound_CMx_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMy_Inv() const final { return AllBound_CMy_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMz_Inv() const final { return AllBound_CMz_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFx_Inv() const final { return AllBound_CFx_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFy_Inv() const final { return AllBound_CFy_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFz_Inv() const final { return AllBound_CFz_Inv; }

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
   * \param[in]  config       - Definition of the particular problem.
   * \param[in]  nFaceSimul   - Number of faces that are treated simultaneously
                                to improve performance.
   * \param[in]  NPad         - Value of the padding parameter to obtain optimal
                                performance in the gemm computations.
   * \param[in]  nPoints      - Number of points per face for which the fluxes
                                must be computed.
   * \param[in]  normalsFace  - The normals in the points for the faces.
   * \param[in]  gridVelsFace - The grid velocities in the points for the faces.
   * \param[in]  solL         - Solution in the left state of the points.
   * \param[in]  solR         - Solution in the right state of the points.
   * \param[out] fluxes       - Inviscid fluxes in the points.
   * \param[in]  numerics     - Object, which contains the Riemann solver.
   */
  void ComputeInviscidFluxesFace(CConfig              *config,
                                 const unsigned short nFaceSimul,
                                 const unsigned short NPad,
                                 const unsigned long  nPoints,
                                 const su2double      *normalsFace[],
                                 const su2double      *gridVelsFace[],
                                 const su2double      *solL,
                                 const su2double      *solR,
                                 su2double            *fluxes,
                                 CNumerics            *numerics);

  /*!
   * \brief Function, which computes the inviscid fluxes in the face integration
            points of a chunk of matching internal faces.
   * \param[in]  config   - Definition of the particular problem.
   * \param[in]  lBeg     - Start index in matchingInternalFaces for which
                            the inviscid fluxes should be computed.
   * \param[in]  lEnd     - End index (not included) in matchingInternalFaces
                            for which the inviscid fluxes should be computed.
   * \param[in]  NPad     - Value of the padding parameter to obtain optimal
                            performance in the gemm computations.
   * \param[out] solIntL  - Solution in the left state of the integration points.
   * \param[out] solIntR  - Solution in the right state of the integration points.
   * \param[out] fluxes   - Inviscid fluxes in the integration points.
   * \param[in]  numerics - Object, which contains the Riemann solver.
   */
  void InviscidFluxesInternalMatchingFace(CConfig              *config,
                                          const unsigned long  lBeg,
                                          const unsigned long  lEnd,
                                          const unsigned short NPad,
                                          su2double            *solIntL,
                                          su2double            *solIntR,
                                          su2double            *fluxes,
                                          CNumerics            *numerics);
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
   * \brief Function, which computes the boundary states in the integration points
            of the boundary face by applying the inviscid wall boundary conditions.
   * \param[in]  config     - Definition of the particular problem.
   * \param[in]  nFaceSimul - Number of faces that are treated simultaneously
                              to improve performance.
   * \param[in]  NPad       - Value of the padding parameter to obtain optimal
                              performance in the gemm computations.
   * \param[in]  surfElem   - Surface boundary elements for which the left state must
                              be computed.
   * \param[in]  solIntL    - Left states in the integration points of the face.
   * \param[out] solIntR    - Right states in the integration points of the face.
   */
  void BoundaryStates_Euler_Wall(CConfig                  *config,
                                 const unsigned short     nFaceSimul,
                                 const unsigned short     NPad,
                                 const CSurfaceElementFEM *surfElem,
                                 const su2double          *solIntL,
                                 su2double                *solIntR);

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
                                                   CVolumeElementFEM    *elem,
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
                                                   CVolumeElementFEM    *elem,
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
                                                      CVolumeElementFEM    *elem,
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
                                                      CVolumeElementFEM    *elem,
                                                      const su2double      *sol,
                                                      const unsigned short nSimul,
                                                      const unsigned short NPad,
                                                      su2double            *res,
                                                      su2double            *work);

  /*!
   * \brief Function, which computes the graph of the spatial discretization
            for the locally owned DOFs.
   * \param[in] DGGeometry - Geometrical definition of the DG problem.
   * \param[in] config     - Definition of the particular problem.
   */
  void DetermineGraphDOFs(const CMeshFEM *FEMGeometry,
                          CConfig        *config);

  /*!
   * \brief Function, which determines the meta data needed for the computation
            of the Jacobian of the spatial residual.
   * \param[in] DGGeometry     - Geometrical definition of the DG problem.
   * \param[in] colorLocalDOFs - Color of the locally stored DOFs.
   */
  void MetaDataJacobianComputation(const CMeshFEM    *FEMGeometry,
                                   const vector<int> &colorLocalDOFs);

  /*!
   * \brief Function, which sets up the list of tasks to be carried out in the
            computationally expensive part of the solver.
   * \param[in] config - Definition of the particular problem.
   */
  void SetUpTaskList(CConfig *config);

  /*!
   * \brief Function, which sets up the persistent communication of the flow
            variables in the DOFs.
   * \param[in] DGGeometry - Geometrical definition of the DG problem.
   * \param[in] config     - Definition of the particular problem.
   */
  void Prepare_MPI_Communication(const CMeshFEM *FEMGeometry,
                                 CConfig *config);

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
   * \param[in]  config    - Definition of the particular problem.
   * \param[in]  useADER   - Whether or not the ADER residual must be multiplied.
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   * \param[out] workArray - Work array.
   */
  void MultiplyResidualByInverseMassMatrix(CConfig             *config,
                                           const bool          useADER,
                                           const unsigned long elemBeg,
                                           const unsigned long elemEnd,
                                           su2double           *workArray);

  /*!
   * \brief Function, which computes the residual contribution from a boundary
            face in an inviscid computation when the boundary conditions have
            already been applied.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     nFaceSimul    - Number of faces that are treated simultaneously
                                    to improve performance.
   * \param[in]     NPad          - Value of the padding parameter to obtain optimal
                                    performance in the gemm computations.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     surfElem      - Surface boundary element for which the
                                    contribution to the residual must be computed.
   * \param[in]     solInt0       - Solution in the integration points of side 0.
                                    It is not const, because the array is used for
                                    temporary storage for the residual.
   * \param[in]     solInt1       - Solution in the integration points of side 1.
   * \param[out]    fluxes        - Temporary storage for the fluxes in the
                                    integration points.
   * \param[out]    resFaces      - Array to store the residuals of the face.
   * \param[in,out] indResFaces   - Index in resFaces, where the current residual
                                    should be stored. It is updated in the function
                                    for the next boundary element.
   */
  void ResidualInviscidBoundaryFace(CConfig                  *config,
                                    const unsigned short     nFaceSimul,
                                    const unsigned short     NPad,
                                    CNumerics                *conv_numerics,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *solInt0,
                                    su2double                *solInt1,
                                    su2double                *fluxes,
                                    su2double                *resFaces,
                                    unsigned long            &indResFaces);

protected:
  /*!
   * \brief Template function, which determines some meta data for the chunk of
            elements/faces that must be treated simulaneously.
   * \param[in]  elem       - Const pointer the volume or face elements for which
                              the meta data must be computed.
   * \param[in]  l          - Start index for the current chunk of elements/faces.
   * \param[in]  elemEnd    - End index (index not included) of the elements to be
                              treated in the residual computation from which this
                              function is called.
   * \param[in]  nElemSimul - Desired number of elements/faces that must be treated
                              simultaneously for optimal performance.
   * \param[in]  nPadMin    - Minimum number of the padding value in the gemm calls.
   * \param[out] lEnd       - Actual end index (not included) for this chunk of
                              elements.
   * \param[out] ind        - Index in the standard elements to which this chunk of
                              elements can be mapped.
   * \param[out] llEnd      - Actual number of elements/faces that are treated
                              simultaneously, llEnd = lEnd - l.
   * \param[out] NPad       - Actual padded N value in the gemm computations for
                              this chunk of elements.
   */
  template <class TElemType>
  void MetaDataChunkOfElem(const TElemType      *elem,
                           const unsigned long  l,
                           const unsigned long  elemEnd,
                           const unsigned short nElemSimul,
                           const unsigned short nPadMin,
                           unsigned long        &lEnd,
                           unsigned short       &ind,
                           unsigned short       &llEnd,
                           unsigned short       &NPad) {

    /* Determine the end index for this chunk of elements that must be
       treated simulaneously. The elements of this chunk must have the
       same standard element in order to make this work. */
    const unsigned long lEndMax = min(l+nElemSimul, elemEnd);

    ind = elem[l].indStandardElement;
    for(lEnd=l+1; lEnd<lEndMax; ++lEnd) {
      if(elem[lEnd].indStandardElement != ind) break;
    }

    /* Store the number of elements that are treated simultaneously in this chunk
       in llEnd and determine the padded N value in the gemm computations. */
    llEnd = lEnd - l;
    NPad  = llEnd*nVar;
    if( NPad%nPadMin ) NPad += nPadMin - (NPad%nPadMin);
  }
};
