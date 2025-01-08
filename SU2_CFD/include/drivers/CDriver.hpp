/*!
 * \file CDriver.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
 * \author T. Economon, H. Kline, R. Sanchez
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

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/parallelization/mpi_structure.hpp"
#include "../integration/CIntegration.hpp"
#include "../interfaces/CInterface.hpp"
#include "../solvers/CSolver.hpp"
#include "CDriverBase.hpp"

using namespace std;

class CInterpolator;
class CIteration;
class COutput;

/*!
 * \class CDriver
 * \ingroup Drivers
 * \brief Parent class for driving an iteration of a single or multi-zone problem.
 * \author T. Economon
 */

class CDriver : public CDriverBase {
 protected:
  su2double
      UsedTimeOutput; /*!< \brief Elapsed time between Start and Stop point of the timer for tracking output phase.*/

  su2double BandwidthSum =
      0.0;                    /*!< \brief Aggregate value of the bandwidth for writing restarts (to be average later).*/
  unsigned long IterCount,    /*!< \brief Iteration count stored for performance benchmarking.*/
      OutputCount;            /*!< \brief Output count stored for performance benchmarking.*/
  unsigned long DOFsPerPoint; /*!< \brief Number of unknowns at each vertex, i.e., number of equations solved. */
  su2double Mpoints; /*!< \brief Total number of grid points in millions in the calculation (including ghost points).*/
  su2double
      MpointsDomain; /*!< \brief Total number of grid points in millions in the calculation (excluding ghost points).*/
  su2double MDOFs;   /*!< \brief Total number of DOFs in millions in the calculation (including ghost points).*/
  su2double MDOFsDomain; /*!< \brief Total number of DOFs in millions in the calculation (excluding ghost points).*/

  bool StopCalc,   /*!< \brief Stop computation flag.*/
      mixingplane, /*!< \brief mixing-plane simulation flag.*/
      fsi,         /*!< \brief FSI simulation flag.*/
      fem_solver;  /*!< \brief FEM fluid solver simulation flag. */

  CFreeFormDefBox*** FFDBox;            /*!< \brief FFD FFDBoxes of the problem. */

  CIteration*** iteration_container;      /*!< \brief Container vector with all the iteration methods. */
  CIntegration**** integration_container; /*!< \brief Container vector with all the integration methods. */
  vector<vector<unique_ptr<CInterpolator>>>
      interpolator_container; /*!< \brief Definition of the interpolation method between non-matching discretizations of
                                 the interface. */
  CInterface*** interface_container; /*!< \brief Definition of the interface of information and physics. */
  bool dry_run;                      /*!< \brief Flag if SU2_CFD was started as dry-run via "SU2_CFD -d <config>.cfg" */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   * \param[in] dummy_geo - Dummy geometric definition of the problem.
   */
  CDriver(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator, bool dummy_geo);

  /*!
   * \brief Destructor of the class.
   */
  ~CDriver(void) override;

  /*!
   * \brief A virtual member.
   */
  void Run() override{};

 protected:
  /*!
   * \brief Initialize containers.
   */
  void InitializeContainers();

  /*!
   * \brief Read in the config and mesh files.
   * \param[in] config - Definition of the particular problem.
   * \param[in] driver_config - Definition of the driver configuration.
   */
  void PreprocessInput(CConfig**& config, CConfig*& driver_config);

  /*!
   * \brief Construction of the edge-based data structure and the multi-grid structure.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] dummy - Definition of the dummy driver.
   */
  void InitializeGeometry(CConfig* config, CGeometry**& geometry, bool dummy);

  /*!
   * \brief Do the geometrical preprocessing for the DG FEM solver.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void InitializeGeometryDGFEM(CConfig* config, CGeometry**& geometry);

  /*!
   * \brief InitializeGeometryFVM
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void InitializeGeometryFVM(CConfig* config, CGeometry**& geometry);

  /*!
   * \brief Definition of the physics iteration class or within a single zone.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iteration - Pointer to the iteration container to be instantiated.
   */
  void PreprocessIteration(CConfig* config, CIteration*& iteration) const;

  /*!
   * \brief Definition and allocation of all solution classes.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   */
  void InitializeSolver(CConfig* config, CGeometry** geometry, CSolver***& solver);

  /*!
   * \brief Preprocess the inlets via file input for all solvers.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void PreprocessInlet(CSolver*** solver, CGeometry** geometry, CConfig* config) const;

  /*!
   * \brief Restart of the solvers from the restart files.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] update_geo - Boolean to indicate if geometry should be updated.
   */
  void RestartSolver(CSolver*** solver, CGeometry** geometry, CConfig* config, bool update_geo);

  /*!
   * \brief Definition and allocation of all solution classes.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iInst - Current solver instance.
   */
  void FinalizeSolver(CSolver**** solver, CGeometry** geometry, CConfig* config, unsigned short val_iInst);

  /*!
   * \brief Definition and allocation of all integration classes.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] integration - Container vector with all the integration methods.
   */
  void InitializeIntegration(CConfig* config, CSolver** solver, CIntegration**& integration) const;

  /*!
   * \brief Definition and allocation of all integration classes.
   * \param[in] integration - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iInst - Current solver instance.
   */
  void FinalizeIntegration(CIntegration*** integration, CGeometry** geometry, CConfig* config,
                                  unsigned short val_iInst);

  /*!
   * \brief Definition and allocation of all interface classes.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] interface_types - Type of coupling between the distinct (physical) zones.
   * \param[in] interface - Class defining the physical transfer of information.
   * \param[in] interpolation -  Object defining the interpolation.
   */
  void InitializeInterface(CConfig** config, CSolver***** solver, CGeometry**** geometry,
                               unsigned short** interface_types, CInterface*** interface,
                               vector<vector<unique_ptr<CInterpolator>>>& interpolation);

  /*!
   * \brief Definition and allocation of all solver classes.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   */
  void InitializeNumerics(CConfig* config, CGeometry** geometry, CSolver*** solver, CNumerics****& numerics) const;

  /*!
   * \brief Helper to instantiate turbulence numerics specialized for different flow solvers.
   */
  template <class FlowIndices>
  void InstantiateTurbulentNumerics(unsigned short nVar_Turb, int offset, const CConfig* config,
                                    const CSolver* turb_solver, CNumerics****& numerics) const;

  /*!
   * \brief Helper to instantiate transition numerics specialized for different flow solvers.
   */
  template <class FlowIndices>
  void InstantiateTransitionNumerics(unsigned short nVar_Trans, int offset, const CConfig* config,
                                     const CSolver* trans_solver, CNumerics****& numerics) const;
  /*!
   * \brief Helper to instantiate species transport numerics specialized for different flow solvers.
   */
  template <class FlowIndices>
  void InstantiateSpeciesNumerics(unsigned short nVar_Species, int offset, const CConfig* config,
                                  const CSolver* species_solver, CNumerics****& numerics) const;

  /*!
   * \brief Definition and allocation of all solver classes.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iInst - Current solver instance.
   */
  void FinalizeNumerics(CNumerics***** numerics, CSolver*** solver, CGeometry** geometry, CConfig* config,
                               unsigned short val_iInst);

  /*!
   * \brief GridMovement_Preprocessing
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] iteration - Container vector with all the iteration methods.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   */
  void PreprocessDynamicMesh(CConfig* config, CGeometry** geometry, CSolver*** solver, CIteration* iteration,
                                 CVolumetricMovement*& grid_movement, CSurfaceMovement*& surface_movement) const;

  /*!
   * \brief Initialize Python interface functionalities. When using multigrid,
   * it is important to call this after modifying custom boundary values.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   */
  void PreprocessPythonInterface(CConfig** config, CGeometry**** geometry, CSolver***** solver);

  /*!
   * \brief Preprocess the output container.
   * \param[in] config - Definition of the particular problem.
   * \param[in] driver_config - Definition of the driver configuration.
   * \param[in] output_container - Container vector with all the outputs.
   * \param[in] driver_output - Definition of the driver output.
   */
  void PreprocessOutput(CConfig** config, CConfig* driver_config, COutput**& output_container,
                            COutput*& driver_output);

  /*!
   * \brief Initiate value for static mesh movement such as the gridVel for the ROTATING frame.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void PreprocessStaticMesh(const CConfig* config, CGeometry** geometry);

  /*!
   * \brief Initiate value for static mesh movement such as the gridVel for the ROTATING frame.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] interface - Class defining the physical transfer of information.
   */
  void PreprocessTurbomachinery(CConfig** config, CGeometry**** geometry, CSolver***** solver,
                                    CInterface*** interface);

  /*!
   * \brief Ramp some simulation settings for turbomachinery problems.
   * \param[in] iter - Iteration for the ramp (can be outer or time depending on type of simulation).
   * \note TODO This is not compatible with inner iterations because they are delegated to the iteration class.
   */
  void RampTurbomachineryValues(unsigned long iter);

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone in which the displacements will be predicted.
   * \param[in] targetZone - zone which receives the predicted displacements.
   */
  virtual void PredictDisplacements(unsigned short donorZone, unsigned short targetZone) {}

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone in which the tractions will be predicted.
   * \param[in] targetZone - zone which receives the predicted traction.
   */
  virtual void PredictTractions(unsigned short donorZone, unsigned short targetZone) {}

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone in which the displacements will be transferred.
   * \param[in] targetZone - zone which receives the tractions transferred.
   */
  virtual void TransferDisplacements(unsigned short donorZone, unsigned short targetZone) {}

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone from which the tractions will be transferred.
   * \param[in] targetZone - zone which receives the tractions transferred.
   */
  virtual void TransferTractions(unsigned short donorZone, unsigned short targetZone) {}

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - origin of the information.
   * \param[in] targetZone - destination of the information.
   * \param[in] iOuterIter - Fluid-Structure Interaction subiteration.
   */
  virtual void RelaxationDisplacements(unsigned short donorZone, unsigned short targetZone, unsigned long iOuterIter) {}

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - origin of the information.
   * \param[in] targetZone - destination of the information.
   * \param[in] iOuterIter - Fluid-Structure Interaction subiteration.
   */
  virtual void RelaxationTractions(unsigned short donorZone, unsigned short targetZone, unsigned long iOuterIter) {}

  /*!
   * \brief Print out the direct residuals.
   * \param[in] kind_recording - Type of recording (full list in ENUM_RECORDING, option_structure.hpp)
   */
  void PrintDirectResidual(RECORDING kind_recording);

  /*!
   * \brief Set the solution of all solvers (adjoint or primal) in a zone.
   * \param[in] iZone - Index of the zone.
   * \param[in] adjoint - True to consider adjoint solvers instead of primal.
   * \param[in] solution - Solution object with interface (iPoint,iVar).
   * \tparam Old - If true set "old solutions" instead.
   */
  template <class Container, bool Old = false>
  void SetAllSolutions(unsigned short iZone, bool adjoint, const Container& solution) {
    const auto nPoint = geometry_container[iZone][INST_0][MESH_0]->GetnPoint();
    for (auto iSol = 0u, offset = 0u; iSol < MAX_SOLS; ++iSol) {
      auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
      if (!(solver && (solver->GetAdjoint() == adjoint))) continue;
      for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint)
        for (auto iVar = 0ul; iVar < solver->GetnVar(); ++iVar)
          if (!Old) {
            solver->GetNodes()->SetSolution(iPoint, iVar, solution(iPoint, offset + iVar));
          } else {
            solver->GetNodes()->SetSolution_Old(iPoint, iVar, solution(iPoint, offset + iVar));
          }
      offset += solver->GetnVar();
    }
  }

  /*!
   * \brief Set the "old solution" of all solvers (adjoint or primal) in a zone.
   */
  template <class Container>
  void SetAllSolutionsOld(unsigned short iZone, bool adjoint, const Container& solution) {
    SetAllSolutions<Container, true>(iZone, adjoint, solution);
  }

  /*!
   * \brief Get the solution of all solvers (adjoint or primal) in a zone.
   * \param[in] iZone - Index of the zone.
   * \param[in] adjoint - True to consider adjoint solvers instead of primal.
   * \param[out] solution - Solution object with interface (iPoint,iVar).
   */
  template <class Container>
  void GetAllSolutions(unsigned short iZone, bool adjoint, Container& solution) const {
    const auto nPoint = geometry_container[iZone][INST_0][MESH_0]->GetnPoint();
    for (auto iSol = 0u, offset = 0u; iSol < MAX_SOLS; ++iSol) {
      auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
      if (!(solver && (solver->GetAdjoint() == adjoint))) continue;
      const auto& sol = solver->GetNodes()->GetSolution();
      for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint)
        for (auto iVar = 0ul; iVar < solver->GetnVar(); ++iVar)
          solution(iPoint, offset + iVar) = SU2_TYPE::GetValue(sol(iPoint, iVar));
      offset += solver->GetnVar();
    }
  }

  /*!
   * \brief Sum the number of primal or adjoint variables for all solvers in a given zone.
   * \param[in] iZone - Index of the zone.
   * \param[in] adjoint - True to consider adjoint solvers instead of primal.
   * \return Total number of solution variables.
   */
  unsigned short GetTotalNumberOfVariables(unsigned short iZone, bool adjoint) const {
    unsigned short nVar = 0;
    for (auto iSol = 0u; iSol < MAX_SOLS; iSol++) {
      auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
      if (solver && (solver->GetAdjoint() == adjoint)) nVar += solver->GetnVar();
    }
    return nVar;
  }

 public:
  /*!
   * \brief Launch the computation for all zones and all physics.
   */
  virtual void StartSolver() {}

  /*!
   * \brief Deallocation routine
   */
  void Finalize() override;

  /*!
   * \brief Perform some pre-processing before an iteration of the physics.
   */
  virtual void Preprocess(unsigned long TimeIter) {}

  /*!
   * \brief Monitor the computation.
   */
  virtual bool Monitor(unsigned long TimeIter) { return false; }

  /*!
   * \brief Output the solution in solution file.
   */
  virtual void Output(unsigned long TimeIter) {}

  /*!
   * \brief Perform a dynamic mesh deformation, including grid velocity computation and update of the multi-grid
   * structure.
   */
  virtual void DynamicMeshUpdate(unsigned long TimeIter) {}

  /*!
   * \brief Update the dual-time solution.
   */
  virtual void Update() {}

  /*!
   * \brief Perform a mesh deformation as initial condition.
   */
  virtual void SetInitialMesh() {}

/// \addtogroup PySU2
/// \{

  /*!
   * \brief Process the boundary conditions and update the multi-grid structure.
   */
  void BoundaryConditionsUpdate();

  /*!
   * \brief Get the number of time iterations.
   * \return Number of time iterations.
   */
  unsigned long GetNumberTimeIter() const;

  /*!
   * \brief Get the current time iteration.
   * \return Current time iteration.
   */
  unsigned long GetTimeIter() const;

  /*!
   * \brief Get the unsteady time step.
   * \return Unsteady time step.
   */
  passivedouble GetUnsteadyTimeStep() const;

  /*!
   * \brief Get the name of the output file for the surface.
   * \return File name for the surface output.
   */
  string GetSurfaceFileName() const;

  /*!
   * \brief Set the position of the heat source.
   * \param[in] alpha - Angle of rotation respect to Z axis.
   * \param[in] pos_x - Position X.
   * \param[in] pos_y - Position Y.
   * \param[in] pos_z - Position Z.
   */
  void SetHeatSourcePosition(passivedouble alpha, passivedouble pos_x, passivedouble pos_y, passivedouble pos_z);

  /*!
   * \brief Set the direction of the inlet.
   * \param[in] iMarker - Marker index.
   * \param[in] alpha - Angle (Zpos).
   */
  void SetInletAngle(unsigned short iMarker, passivedouble alpha);

  /*!
   * \brief Set the angle of attack of the farfield.
   * \param[in] alpha - Angle (degree).
   */
  void SetFarFieldAoA(passivedouble alpha);

  /*!
   * \brief Set the angle of sideslip of the farfield.
   * \param[in] beta - Angle (degree).
   */
  void SetFarFieldAoS(passivedouble beta);

  /*!
   * \brief Set the dynamic mesh translation rates.
   * \param[in] xDot - Value of translational velocity in x-direction.
   * \param[in] yDot - Value of translational velocity in y-direction.
   * \param[in] zDot - Value of translational velocity in z-direction.
   */
  void SetTranslationRate(passivedouble xDot, passivedouble yDot, passivedouble zDot);

  /*!
   * \brief Set the dynamic mesh rotation rates.
   * \param[in] rot_x - Value of Angular velocity about x-axes.
   * \param[in] rot_y - Value of Angular velocity about y-axes.
   * \param[in] rot_z - Value of Angular velocity about z-axes.
   */
  void SetRotationRate(passivedouble rot_x, passivedouble rot_y, passivedouble rot_z);

/// \}
};

/*!
 * \class CFluidDriver
 * \ingroup Drivers
 * \brief Class for driving an iteration of the physics within multiple zones.
 * \author T. Economon, G. Gori
 */
class CFluidDriver : public CDriver {
 protected:
  unsigned long Max_Iter;

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CFluidDriver(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator);

  /*!
   * \brief Transfer data among different zones (multiple zone).
   */
  void TransferData(unsigned short donorZone, unsigned short targetZone);

 public:
  /*!
   * \brief Destructor of the class.
   */
  ~CFluidDriver(void) override;

  /*!
   * \brief Launch the computation for all zones and all physics.
   */
  void StartSolver() override;

  /*!
   * \brief Run a single iteration of the physics within multiple zones.
   */
  void Run() override;

  /*!
   * \brief Update the dual-time solution within multiple zones.
   */
  void Update() override;

  /*!
   * \brief Output the solution in solution file.
   */
  void Output(unsigned long InnerIter) override;

  /*!
   * \brief Monitor the computation.
   */
  bool Monitor(unsigned long ExtIter) override;

  /*!
   * \brief Perform some pre-processing before an iteration of the physics.
   */
  void Preprocess(unsigned long Iter) override;

  /*!
   * \brief Perform a dynamic mesh deformation, included grid velocity computation and the update of the multi-grid
   * structure (multiple zone).
   */
  void DynamicMeshUpdate(unsigned long TimeIter) override;
};

/*!
 * \class CHBDriver
 * \ingroup Drivers
 * \brief Class for driving an iteration of Harmonic Balance (HB) method problem using multiple time zones.
 * \author T. Economon
 */
class CHBDriver : public CFluidDriver {
 private:
  unsigned short nInstHB;
  su2double** D; /*!< \brief Harmonic Balance operator. */

  /*!
   * \brief Computation and storage of the Harmonic Balance method source terms.
   * \author T. Economon, K. Naik
   * \param[in] iZone - Current zone number.
   */
  void SetHarmonicBalance(unsigned short iZone);

  /*!
   * \brief Precondition Harmonic Balance source term for stability
   * \author J. Howison
   */
  void StabilizeHarmonicBalance();

  /*!
   * \brief Computation of the Harmonic Balance operator matrix for harmonic balance.
   * \author A. Rubino, S. Nimmagadda
   */
  void ComputeHBOperator();

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CHBDriver(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CHBDriver(void) override;

  /*!
   * \brief Run a single iteration of a Harmonic Balance problem.
   */
  void Run() override;

  /*!
   * \brief Update the solution for the Harmonic Balance.
   */
  void Update() override;
};
