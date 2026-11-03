/*!
 * \file CTurbSASolver.hpp
 * \brief Headers of the CTurbSASolver class
 * \author A. Bueno.
 * \version 8.6.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "CTurbSolver.hpp"

/*!
 * \class CTurbSASolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */

class CTurbSASolver final : public CTurbSolver {
private:

  su2double nu_tilde_Engine[4] = {0.0};
  su2double nu_tilde_ActDisk[4] = {0.0};

  /*--- FIML Neural Network Components ---*/
  su2double*** weights = nullptr;              /*!< \brief NN weights [layer][neuron_from][neuron_to] */
  su2double*** weight_gradients = nullptr;     /*!< \brief Weight gradients [layer][neuron_from][neuron_to] */
  unsigned short n_neurons = 0;                /*!< \brief Number of neurons per hidden layer */
  unsigned short n_hidden_layers = 0;          /*!< \brief Number of hidden layers */
  unsigned short n_inputs = 4;                 /*!< \brief Number of NN inputs (features) */

  su2double box_cox_lambda[4] = {1.0, 1.0, 1.0, 1.0};  /*!< \brief Box-Cox transform parameters */
  su2double feature_mean[4] = {0.0, 0.0, 0.0, 0.0};    /*!< \brief Feature means for scaling */
  su2double feature_std[4] = {1.0, 1.0, 1.0, 1.0};     /*!< \brief Feature std devs for scaling */
  const su2double feature_shift[4] = {0.0, 0.0, 0.5, 0.0};  /*!< \brief Shift that makes each feature non-negative before Box-Cox. */
  bool nn_scaling_set = false;                 /*!< \brief Scaling computed or loaded with the weights. */
  bool targets_loaded = false;                 /*!< \brief Beta targets have been read. */
  unsigned short epoch_counter = 0;            /*!< \brief Completed training epochs. */

  bool filter_shield = false;                  /*!< \brief Apply spatial filtering */
  su2double Total_CpDiff_FIML = 0.0;          /*!< \brief FIML pressure coefficient objective */
  su2double Total_ClDiff_FIML = 0.0;          /*!< \brief FIML lift coefficient objective */
  su2double Total_CdDiff_FIML = 0.0;          /*!< \brief FIML drag coefficient objective */

  /*!
   * \brief Initialize neural network weights and structure.
   * \param[in] config - Definition of the particular problem.
   */
  void InitializeNeuralNetwork(const CConfig* config);

  /*!
   * \brief Forward propagate through neural network to compute beta_fiml at all points.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver_container - Container with all solvers.
   * \param[in] geometry - Geometrical definition.
   */
  void ForwardPropagate(const CConfig* config, CSolver** solver_container, CGeometry* geometry);

  /*!
   * \brief Compute neural network input features at a point.
   * \param[in] iPoint - Grid point index.
   * \param[in] solver_container - Container with all solvers.
   * \param[in] geometry - Geometrical definition.
   * \param[out] features - Array of 4 feature values.
   */
  void ComputeNNInputFeatures(unsigned long iPoint, CSolver** solver_container,
                               CGeometry* geometry, su2double* features) const;

  /*!
   * \brief Apply Box-Cox transform and z-score normalization to input features.
   * \param[in,out] features - Feature array to scale (in place).
   */
  void ScaleNNInputs(su2double* features) const;

  /*!
   * \brief Check if point should be filtered based on flow physics.
   * \param[in] iPoint - Grid point index.
   * \return True if point should be filtered out.
   */
  bool ApplyFilterShield(unsigned long iPoint) const;

  /*!
   * \brief Fit Box-Cox lambda per feature by maximum likelihood, and the mean and std of the transformed features.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver_container - Container with all solvers.
   * \param[in] geometry - Geometrical definition.
   */
  void ComputeBoxCoxLambda(const CConfig* config, CSolver** solver_container, CGeometry* geometry);

  /*!
   * \brief Backward propagate to compute weight gradients from beta targets.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver_container - Container with all solvers.
   * \param[in] geometry - Geometrical definition.
   */
  void BackwardPropagate(const CConfig* config, CSolver** solver_container, CGeometry* geometry);

  /*!
   * \brief Update neural network weights using gradient descent.
   * \param[in] learning_rate - Learning rate for weight update.
   */
  void UpdateWeights(su2double learning_rate);

  /*!
   * \brief Load beta target values for training.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition.
   */
  void LoadBetaTargets(const CConfig* config, CGeometry* geometry);

  /*!
   * \brief Compute training loss (MSE between predicted and target beta).
   * \param[in] geometry - Geometrical definition.
   * \return Total training loss.
   */
  su2double ComputeTrainingLoss(CGeometry* geometry) const;

  /*!
   * \brief Save neural network weights to file.
   * \param[in] filename - Output filename for weights.
   */
  void SaveNeuralNetworkWeights(const string& filename = "nn_weights.dat") const;

  /*!
   * \brief Load neural network weights from file.
   * \param[in] filename - Input filename containing weights.
   * \return True if weights were successfully loaded, false otherwise.
   */
  bool LoadNeuralNetworkWeights(const string& filename = "nn_weights.dat");

  /*!
   * \brief Load beta_fiml from the restart field Beta_FIML when present.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void LoadRestartAuxiliaryFields(const CGeometry* geometry, const CConfig* config) override;

  /*!
   * \brief A virtual member.
   * \param[in] solver - Solver container
   * \param[in] geometry - Geometrical definition.
   * \param[in] config - Definition of the particular problem.
   */
  void SetDES_LengthScale(CSolver** solver,
                          CGeometry *geometry,
                          CConfig *config);

  /*!
   * \brief Mark the points that are located inside the box where the Stochastic Backscatter Model is active.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition.
   */
  void SetBackscatterInBox(CConfig *config, CGeometry* geometry);

  /*!
   * \brief Update the source terms of the stochastic equations (Stochastic Backscatter Model).
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition.
   */
  void SetLangevinSourceTerms(CConfig *config, CGeometry* geometry);

  /*!
   * \brief Apply Laplacian smoothing to the source terms in Langevin equations (Stochastic Backscatter Model).
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition.
   */
  void SmoothLangevinSourceTerms(CConfig* config, CGeometry* geometry);

  /*!
   * \brief Compute nu tilde from the wall functions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void SetTurbVars_WF(CGeometry *geometry,
                     CSolver **solver_container,
                     const CConfig *config,
                     unsigned short val_marker);

  /*!
   * \brief Compute a suitable under-relaxation parameter to limit the change in the solution variables over
   * a nonlinear iteration for stability.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeUnderRelaxationFactor(CSolver** solver_container, const CConfig *config) final;

  /*!
   * \brief Resolve the compile-time parameters of CScalarFlux_SA and run one of this solver's
   *        boundaries through the shared boundary flux pass.
   * \param[in] opt - Flags of the boundary, from one of ScalarFluxOptions' named constructors.
   */
  void BoundaryFlux(CGeometry* geometry, CSolver** solver_container, CConfig* config, const ScalarFluxOptions& opt,
                    unsigned short val_marker);

public:
  /*!
   * \brief Constructor.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] FluidModel
   */
  CTurbSASolver(CGeometry *geometry, CConfig *config, const CSolver* flow_solver, unsigned short iMesh,
                CFluidModel* FluidModel);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurbSASolver() override;

  /*!
   * \brief Restart residual and compute gradients.
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
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Postprocessing(CGeometry *geometry,
                      CSolver **solver_container,
                      CConfig *config,
                      unsigned short iMesh) override;

  /*!
   * \brief Compute the spatial integration using the CScalarFlux_SA edge kernel, which computes
   *        and writes both the convective and the diffusive term of every edge.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Unused, kept only for the boundary conditions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Upwind_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics** numerics_container,
                       CConfig* config, unsigned short iMesh) override;

  /*!
   * \brief Impose the Far Field boundary condition, via the CScalarFlux_SA edge kernel.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Unused, kept only for the boundary condition dispatch.
   * \param[in] visc_numerics - Unused, kept only for the boundary condition dispatch.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) override;

  /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Residual(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics **numerics_container,
                       CConfig *config,
                       unsigned short iMesh) override;

  /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Template(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics *numerics,
                       CConfig *config,
                       unsigned short iMesh) override;

  /*!
   * \brief Impose the Navier-Stokes wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatFlux_Wall(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) override;

  /*!
   * \brief Impose the Navier-Stokes wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Isothermal_Wall(CGeometry *geometry,
                          CSolver **solver_container,
                          CNumerics *conv_numerics,
                          CNumerics *visc_numerics,
                          CConfig *config,
                          unsigned short val_marker) override;

  /*!
   * \brief Impose the inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet(CGeometry *geometry,
                CSolver **solver_container,
                CNumerics *conv_numerics,
                CNumerics *visc_numerics,
                CConfig *config,
                unsigned short val_marker) override;

  /*!
   * \brief Impose the inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet_Turbo(CGeometry *geometry,
                      CSolver **solver_container,
                      CNumerics *conv_numerics,
                      CNumerics *visc_numerics,
                      CConfig *config,
                      unsigned short val_marker) override;

  /*!
   * \brief Impose the inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet_MixingPlane(CGeometry *geometry,
                            CSolver **solver_container,
                            CNumerics *conv_numerics,
                            CNumerics *visc_numerics,
                            CConfig *config,
                            unsigned short val_marker) override;

  /*!
   * \brief Impose the outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Outlet(CGeometry *geometry,
                 CSolver **solver_container,
                 CNumerics *conv_numerics,
                 CNumerics *visc_numerics,
                 CConfig *config,
                 unsigned short val_marker) override;

  /*!
   * \brief Impose the engine inflow boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Engine_Inflow(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) override;

  /*!
   * \brief Impose the engine exhaust boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Engine_Exhaust(CGeometry *geometry,
                         CSolver **solver_container,
                         CNumerics *conv_numerics,
                         CNumerics *visc_numerics,
                         CConfig *config,
                         unsigned short val_marker) override;

  /*!
   * \brief Impose an actuator disk inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ActDisk_Inlet(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) override;

  /*!
   * \brief Impose an actuator disk outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ActDisk_Outlet(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) override;

  /*!
   * \brief Impose an actuator disk inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   * \param[in] val_inlet_surface - Boolean for whether val_marker is an inlet
   */
  void BC_ActDisk(CGeometry *geometry,
                  CSolver **solver_container,
                  CNumerics *conv_numerics,
                  CNumerics *visc_numerics,
                  CConfig *config,
                  unsigned short val_marker,
                  bool val_inlet_surface) override;

  /*!
   * \brief Impose the fluid interface (sliding mesh) boundary condition, via the CScalarFlux_SA
   *        edge kernel. The convective term is a per-donor weighted average, computed in the same
   *        pass that fills the ghost row of each donor; the diffusive term is computed once per
   *        vertex, after the donor loop, from the ghost state the last donor left behind.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Unused, kept only for the boundary condition dispatch.
   * \param[in] visc_numerics - Unused, kept only for the boundary condition dispatch.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                          CNumerics *visc_numerics, CConfig *config) override;

  /*!
   * \brief Store of a set of provided inlet profile values at a vertex.
   * \param[in] val_inlet - vector containing the inlet values for the current vertex.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   * \param[in] iVertex - Vertex of the marker <i>iMarker</i> where the inlet is being set.
   */
  void SetInletAtVertex(const su2double *val_inlet,
                        unsigned short iMarker,
                        unsigned long iVertex) override;

  /*!
   * \brief Get the set of values imposed at an inlet.
   * \param[in] iMarker - Index of the surface marker.
   * \param[in] iVertex - Vertex of the marker <i>iMarker</i> where the inlet is being set.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in,out] val_inlet - vector returning the inlet values for the current vertex.
   * \return Value of the face area at the vertex.
   */
  su2double GetInletAtVertex(unsigned short iMarker, unsigned long iVertex,
                             const CGeometry* geometry, su2double* val_inlet) const override;

  /*!
   * \brief Set a uniform inlet profile
   *
   * The values at the inlet are set to match the values specified for
   * inlets in the configuration file.
   *
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   */
  void SetUniformInlet(const CConfig* config, unsigned short iMarker) override;

  /*!
   * \brief Get the value of nu tilde at the far-field.
   * \return Value of nu tilde at the far-field.
   */
  inline su2double GetNuTilde_Inf(void) const override { return Solution_Inf[0]; }

};
