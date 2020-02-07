/*!
 * \file CTurbSolver.hpp
 * \brief Headers of the CTurbSolver class
 * \author A. Bueno.
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
#include "../variables/CTurbVariable.hpp"

/*!
 * \class CTurbSolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */
class CTurbSolver : public CSolver {
protected:
  su2double *FlowPrimVar_i,    /*!< \brief Store the flow solution at point i. */
  *FlowPrimVar_j,              /*!< \brief Store the flow solution at point j. */
  *lowerlimit,                 /*!< \brief contains lower limits for turbulence variables. */
  *upperlimit;                 /*!< \brief contains upper limits for turbulence variables. */
  su2double Gamma;             /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One;   /*!< \brief Fluids's Gamma - 1.0  . */
  su2double*** Inlet_TurbVars; /*!< \brief Turbulence variables at inlet profiles */

  CTurbVariable* snode;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

  /* Sliding meshes variables */

  su2double ****SlidingState;
  int **SlidingStateNodes;

  CTurbVariable* nodes = nullptr;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() final { return nodes; }

public:

  /*!
   * \brief Constructor of the class.
   */
  CTurbSolver(void);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CTurbSolver(void);

  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTurbSolver(CGeometry* geometry, CConfig *config);

  /*!
   * \brief Compute the spatial integration using a upwind scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */

  void Upwind_Residual(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics *numerics,
                       CConfig *config,
                       unsigned short iMesh) override;

  /*!
   * \brief Compute the viscous residuals for the turbulent equation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void Viscous_Residual(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *numerics,
                        CConfig *config,
                        unsigned short iMesh,
                        unsigned short iRKStep) override;

  /*!
   * \brief Impose the Symmetry Plane boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Sym_Plane(CGeometry      *geometry,
                    CSolver        **solver_container,
                    CNumerics      *conv_numerics,
                    CNumerics      *visc_numerics,
                    CConfig        *config,
                    unsigned short val_marker) override;

  /*!
   * \brief Impose via the residual the Euler wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Euler_Wall(CGeometry      *geometry,
                     CSolver        **solver_container,
                     CNumerics      *conv_numerics,
                     CNumerics      *visc_numerics,
                     CConfig        *config,
                     unsigned short val_marker) override;
  /*!
   * \brief Impose via the residual the Euler wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Riemann(CGeometry *geometry,
                  CSolver **solver_container,
                  CNumerics *conv_numerics,
                  CNumerics *visc_numerics,
                  CConfig *config,
                  unsigned short val_marker) final;

  /*!
   * \brief Impose via the residual the Euler wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_TurboRiemann(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics *conv_numerics,
                       CNumerics *visc_numerics,
                       CConfig *config,
                       unsigned short val_marker) final;

  /*!
   * \brief Impose via the residual the Euler wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Giles(CGeometry *geometry,
                CSolver **solver_container,
                CNumerics *conv_numerics,
                CNumerics *visc_numerics,
                CConfig *config,
                unsigned short val_marker) final;

  /*!
   * \brief Impose a periodic boundary condition by summing contributions from the complete control volume.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Periodic(CGeometry *geometry,
                   CSolver **solver_container,
                   CNumerics *numerics,
                   CConfig *config) final;

  /*!
   * \brief Update the solution using an implicit solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitEuler_Iteration(CGeometry *geometry,
                               CSolver **solver_container,
                               CConfig *config) override;
  /*!
   * \brief Set the total residual adding the term that comes from the Dual Time-Stepping Strategy.
   * \param[in] geometry - Geometric definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void SetResidual_DualTime(CGeometry *geometry,
                            CSolver **solver_container,
                            CConfig *config,
                            unsigned short iRKStep,
                            unsigned short iMesh,
                            unsigned short RunTime_EqSystem) final;

  /*!
   * \brief Compute a suitable under-relaxation parameter to limit the change in the solution variables over a nonlinear iteration for stability.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeUnderRelaxationFactor(CSolver **solver, CConfig *config) final;

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
  * \brief Get the outer state for fluid interface nodes.
  * \param[in] val_marker - marker index
  * \param[in] val_vertex - vertex index
  * \param[in] val_state  - requested state component
  * \param[in] donor_index- index of the donor node to get
  */
  inline su2double GetSlidingState(unsigned short val_marker,
                                   unsigned long val_vertex,
                                   unsigned short val_state,
                                   unsigned long donor_index) const final {
    return SlidingState[val_marker][val_vertex][val_state][donor_index];
  }

  /*!
   * \brief Allocates the final pointer of SlidingState depending on how many donor vertex donate to it. That number is stored in SlidingStateNodes[val_marker][val_vertex].
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   */
  inline void SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex) final {
    int iVar;

    for( iVar = 0; iVar < nVar+1; iVar++){
      if( SlidingState[val_marker][val_vertex][iVar] != NULL )
        delete [] SlidingState[val_marker][val_vertex][iVar];
    }

    for( iVar = 0; iVar < nVar+1; iVar++)
      SlidingState[val_marker][val_vertex][iVar] = new su2double[ GetnSlidingStates(val_marker, val_vertex) ];
  }


  /*!
   * \brief Set the outer state for fluid interface nodes.
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   * \param[in] val_state    - requested state component
   * \param[in] donor_index  - index of the donor node to set
   * \param[in] component    - set value
   */
  inline void SetSlidingState(unsigned short val_marker,
                              unsigned long val_vertex,
                              unsigned short val_state,
                              unsigned long donor_index,
                              su2double component) final {
    SlidingState[val_marker][val_vertex][val_state][donor_index] = component;
  }


  /*!
   * \brief Set the number of outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] value - number of outer states
   */
  inline void SetnSlidingStates(unsigned short val_marker,
                                unsigned long val_vertex,
                                int value) final { SlidingStateNodes[val_marker][val_vertex] = value; }

  /*!
   * \brief Get the number of outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   */
  inline int GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex) const final {
    return SlidingStateNodes[val_marker][val_vertex];
  }

  /*!
   * \brief Set custom turbulence variables at the vertex of an inlet.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] iDim - Index of the turbulence variable (i.e. k is 0 in SST)
   * \param[in] val_turb_var - Value of the turbulence variable to be used.
   */
  inline void SetInlet_TurbVar(unsigned short val_marker,
                               unsigned long val_vertex,
                               unsigned short val_dim,
                               su2double val_turb_var) final {
    /*--- Since this call can be accessed indirectly using python, do some error
     * checking to prevent segmentation faults ---*/
    if (val_marker >= nMarker)
      SU2_MPI::Error("Out-of-bounds marker index used on inlet.", CURRENT_FUNCTION);
    else if (Inlet_TurbVars == NULL || Inlet_TurbVars[val_marker] == NULL)
      SU2_MPI::Error("Tried to set custom inlet BC on an invalid marker.", CURRENT_FUNCTION);
    else if (val_vertex >= nVertex[val_marker])
      SU2_MPI::Error("Out-of-bounds vertex index used on inlet.", CURRENT_FUNCTION);
    else if (val_dim >= nVar)
      SU2_MPI::Error("Out-of-bounds index used for inlet turbulence variable.", CURRENT_FUNCTION);
    else
      Inlet_TurbVars[val_marker][val_vertex][val_dim] = val_turb_var;
  }

};