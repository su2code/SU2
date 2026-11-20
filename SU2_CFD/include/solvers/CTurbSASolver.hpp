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
#include "../fiml/FIMLInputs.hpp"

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

  /*--- FIML neural network. ---*/
  vector<su2activematrix> weights;             /*!< \brief Weights per layer, (inputs + bias) x outputs. */
  vector<su2activematrix> weight_gradients;    /*!< \brief Gradients of the training loss or of the objective. */
  unsigned short n_neurons = 0;                /*!< \brief Number of neurons per hidden layer */
  unsigned short n_hidden_layers = 0;          /*!< \brief Number of hidden layers */
  vector<FIML_INPUT> nn_inputs;                /*!< \brief Network inputs, in order. */
  unsigned short n_inputs = 0;                 /*!< \brief Number of network inputs. */

  su2double box_cox_lambda[FIML::N_INPUTS];    /*!< \brief Box-Cox parameters of the inputs. */
  su2double feature_mean[FIML::N_INPUTS];      /*!< \brief Means of the transformed inputs. */
  su2double feature_std[FIML::N_INPUTS];       /*!< \brief Standard deviations of the transformed inputs. */
  bool nn_scaling_set = false;                 /*!< \brief Scaling computed or loaded with the weights. */
  bool targets_loaded = false;                 /*!< \brief Beta targets have been read. */
  unsigned short epoch_counter = 0;            /*!< \brief Completed training epochs. */

  FIML_SHIELD shield = FIML_SHIELD::NONE;      /*!< \brief Regions where the network is switched off. */
  su2double shield_coeff = 8.0;                /*!< \brief Coefficient of f_d for FIML_SHIELD= FD. */
  su2vector<uint8_t> nn_filtered;              /*!< \brief Points excluded by the filter in the last forward pass. */
  su2double tau_wall_ref = 1.0;                /*!< \brief Flat-plate wall shear stress that normalizes the Delta feature. */
  vector<su2double> adam_m, adam_v;            /*!< \brief Adam moments of the weights, in file order. */
  unsigned long adam_step = 0;                 /*!< \brief Number of Adam updates. */

  /*!
   * \brief Call f(layer, i, j) for every weight, in the order of nn_weights.dat (bias last in each layer).
   */
  template <class F>
  void ForEachWeight(F&& f) const {
    for (unsigned long layer = 0; layer < weights.size(); layer++)
      for (unsigned long i = 0; i < weights[layer].rows(); i++)
        for (unsigned long j = 0; j < weights[layer].cols(); j++) f(layer, i, j);
  }

  /*!
   * \brief Initialize neural network weights and structure.
   * \param[in] config - Definition of the particular problem.
   */
  void InitializeNeuralNetwork(const CConfig* config);

  /*!
   * \brief Activation buffers of one point: each layer's inputs followed by the bias input 1, then tanh(output).
   */
  vector<vector<su2double>> NetworkActivations() const;

  /*!
   * \brief Evaluate the network for scaled inputs, filling the activations of every layer.
   */
  void EvaluateNetwork(const su2double* inputs, vector<vector<su2double>>& activations) const;

  /*!
   * \brief Set beta_fiml at all domain points from the network, or to 1 where the shield excludes the point.
   */
  void ForwardPropagate(const CGeometry* geometry, CSolver** solver_container);

  /*!
   * \brief Local quantities of a point from which the inputs are computed.
   */
  FIML::LocalState GetFIMLLocalState(unsigned long iPoint, const CGeometry* geometry,
                                     const CSolver* const* solver_container) const;

  /*!
   * \brief Unscaled values of the network inputs, in the order of FIML_INPUTS.
   */
  void ComputeNNInputFeatures(const FIML::LocalState& state, su2double* features) const;

  /*!
   * \brief Box-Cox transform and z-score of the inputs that use them; the bounded inputs are left as they are.
   * \param[in,out] features - Input values, scaled in place.
   */
  void ScaleNNInputs(su2double* features) const;

  /*!
   * \brief Whether FIML_SHIELD switches the network off at a point.
   * \return True for LEGACY where DELTA < 0.25 or STRAIN_VORT > 2.
   */
  bool FilterShield(const FIML::LocalState& state) const;

  /*!
   * \brief Set the turbulent flat-plate wall shear stress at half the reference length.
   * \param[in] config - Definition of the particular problem.
   */
  void SetReferenceWallShearStress(const CConfig* config);

  /*!
   * \brief Fit Box-Cox lambda per feature by maximum likelihood, and the mean and std of the transformed features.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver_container - Container with all solvers.
   * \param[in] geometry - Geometrical definition.
   */
  void ComputeBoxCoxLambda(const CConfig* config, CSolver** solver_container, CGeometry* geometry);

  /*!
   * \brief Gradients of the training loss with respect to the weights, averaged over the unfiltered points.
   */
  void BackwardPropagate(const CGeometry* geometry, CSolver** solver_container);

  /*!
   * \brief Update the neural network weights with one Adam step.
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
   * \brief Half the mean square difference between predicted and target beta over the unfiltered points.
   */
  su2double ComputeTrainingLoss() const;

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
   * \brief Objective gradient with respect to the weights, summed over ranks, in file order.
   */
  vector<passivedouble> GlobalNNWeightGradients() const;

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
  ~CTurbSASolver() = default;

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
   * \brief Sum of (beta - 1)^2 over all domain points, the history output BETA_REG.
   * \param[in] config - Definition of the particular problem.
   * \return Global value of the regularization term.
   */
  su2double ComputeFIMLRegularization(const CConfig* config) const override;

  /*!
   * \brief Unscaled values of all FIML inputs at a point, in the order of FIML_INPUT.
   */
  void GetFIMLInputs(unsigned long iPoint, const CGeometry* geometry, const CSolver* const* solver_container,
                     su2double* values) const override;

  /*!
   * \brief Register the neural network weights as AD inputs, or reset their AD indices.
   */
  void RegisterNNWeights(bool reset) override;

  /*!
   * \brief Store the adjoints of the neural network weights.
   */
  void ExtractNNWeightAdjoints() override;

  /*!
   * \brief Write the objective gradient with respect to the weights, summed over ranks, in file order.
   */
  void WriteNNWeightGradients(const string& filename) const override;

  /*!
   * \brief L2 norm of dJ/dbeta (KIND_TRAIN_NN= NONE) or dJ/dw (WEIGHTS), the history output SENS_FIML.
   * \param[in] config - Definition of the particular problem.
   * \return Global norm of the gradient.
   */
  su2double ComputeFIMLSensitivityNorm(const CConfig* config) const override;

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
