/*!
 * \file CTurbSolver.hpp
 * \brief Headers of the CTurbSolver class
 * \author A. Bueno.
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "CScalarSolver.hpp"
#include "../variables/CTurbVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"

/*!
 * \class CTurbSolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */
class CTurbSolver : public CScalarSolver<CTurbVariable> {
protected:

  vector<su2activematrix> Inlet_TurbVars;  /*!< \brief Turbulence variables at inlet profiles */

public:
  /*!
   * \brief Destructor of the class.
   */
  ~CTurbSolver() override;

  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTurbSolver(CGeometry* geometry, CConfig *config, bool conservative);

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
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry** geometry, CSolver*** solver, CConfig* config, int val_iter, bool val_update_geo) override;

  /*!
   * \brief Impose fixed values to turbulence quantities.
   * \details Turbulence quantities are set to far-field values in an upstream half-plane
   * in order to keep them from decaying.
   */
  void Impose_Fixed_Values(const CGeometry *geometry, const CConfig *config) final;

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
    else if (val_vertex >= nVertex[val_marker])
      SU2_MPI::Error("Out-of-bounds vertex index used on inlet.", CURRENT_FUNCTION);
    else if (val_dim >= nVar)
      SU2_MPI::Error("Out-of-bounds index used for inlet turbulence variable.", CURRENT_FUNCTION);
    else
      Inlet_TurbVars[val_marker][val_vertex][val_dim] = val_turb_var;
  }

  /*!
   * \brief Register additional In- or Outputs for RANS.
   * \param[in] input - Boolean whether In- or Output should be registered.
   * \param[in] config - The particular config.
   * \returns The number of extra variables.
   */
  unsigned long RegisterSolutionExtra(bool input, const CConfig* config) final;
  
  /*!
   * \brief Compute the viscous flux for the turbulence equations at a particular edge in a non-conservative manner.
   * \tparam SolverSpecificNumericsTemp - lambda-function, to implement solver specific contributions to numerics.
   * \note The functor has to implement (iPoint, jPoint)
   * \param[in] iEdge - Edge for which we want to compute the flux
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Viscous_Residual(const unsigned long iEdge, const CGeometry* geometry, CSolver** solver_container,
                                     CNumerics* numerics, const CConfig* config) {}
  template<typename SolverSpecificNumericsTemp>
  void TurbViscousResidual(const unsigned long iEdge, const CGeometry* geometry, CSolver** solver_container, 
                                        CNumerics* numerics, const CConfig* config, SolverSpecificNumericsTemp&& SolverSpecificNumerics);

  /*!
   * \brief Compute a suitable under-relaxation parameter to limit the change in the solution variables over
   * a nonlinear iteration for stability.
   * \param[in] allowableRatio - Maximum percentage update in variable per iteration.
   */
  void ComputeUnderRelaxationFactorHelper(su2double allowableRatio);
};

template<typename SolverSpecificNumericsTemp>
void CTurbSolver::TurbViscousResidual(const unsigned long iEdge, const CGeometry* geometry, CSolver** solver_container,
                                      CNumerics* numerics, const CConfig* config, SolverSpecificNumericsTemp&& SolverSpecificNumerics) {
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
    CFlowVariable* flowNodes = solver_container[FLOW_SOL] ?
        su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes()) : nullptr;

  const auto iPoint = geometry->edges->GetNode(iEdge, 0);
  const auto jPoint = geometry->edges->GetNode(iEdge, 1);

  /*--- Lambda function to compute the flux ---*/
  auto ComputeFlux = [&](unsigned long iPoint, unsigned long jPoint, const su2double* normal) {
    numerics->SetCoord(geometry->nodes->GetCoord(iPoint),geometry->nodes->GetCoord(jPoint));
    numerics->SetNormal(normal);

    if (flowNodes) {
      numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), flowNodes->GetPrimitive(jPoint));
    }

  /*--- Solver specific numerics contribution. ---*/
    SolverSpecificNumerics(iPoint, jPoint);

    numerics->SetScalarVar(nodes->GetSolution(iPoint), nodes->GetSolution(jPoint));
    numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(jPoint));

    return numerics->ComputeResidual(config); 
  };

  /*--- Compute fluxes and jacobians i->j ---*/
  const su2double* normal = geometry->edges->GetNormal(iEdge);
  auto residual_ij = ComputeFlux(iPoint, jPoint, normal);

  JacobianScalarType *Block_ii = nullptr, *Block_ij = nullptr, *Block_ji = nullptr, *Block_jj = nullptr;
  if (implicit) {
    Jacobian.GetBlocks(iEdge, iPoint, jPoint, Block_ii, Block_ij, Block_ji, Block_jj);
  }
  if (ReducerStrategy) {
    EdgeFluxes.SubtractBlock(iEdge, residual_ij);
    EdgeFluxesDiff.SetBlock(iEdge, residual_ij);
    if (implicit) {
      /*--- For the reducer strategy the Jacobians are averaged for simplicity. ---*/
      for (int iVar=0; iVar<nVar; iVar++)
        for (int jVar=0; jVar<nVar; jVar++) {
          Block_ij[iVar*nVar + jVar] -= 0.5 * SU2_TYPE::GetValue(residual_ij.jacobian_j[iVar][jVar]);
          Block_ji[iVar*nVar + jVar] += 0.5 * SU2_TYPE::GetValue(residual_ij.jacobian_i[iVar][jVar]);
        }
    }
  } else {
    LinSysRes.SubtractBlock(iPoint, residual_ij);
    if (implicit) {
      for (int iVar=0; iVar<nVar; iVar++)
        for (int jVar=0; jVar<nVar; jVar++) {
          Block_ii[iVar*nVar + jVar] -= SU2_TYPE::GetValue(residual_ij.jacobian_i[iVar][jVar]);
          Block_ij[iVar*nVar + jVar] -= SU2_TYPE::GetValue(residual_ij.jacobian_j[iVar][jVar]);
        }
    }
  }

  /*--- Compute fluxes and jacobians j->i ---*/
  su2double flipped_normal[MAXNDIM];
  for (auto iDim = 0u; iDim < nDim; iDim++) flipped_normal[iDim] = -normal[iDim];

  auto residual_ji = ComputeFlux(jPoint, iPoint, flipped_normal);
  if (ReducerStrategy) {
    EdgeFluxesDiff.AddBlock(iEdge, residual_ji);
    if (implicit) {
      for (int iVar=0; iVar<nVar; iVar++)
        for (int jVar=0; jVar<nVar; jVar++) {
          Block_ij[iVar*nVar + jVar] += 0.5 * SU2_TYPE::GetValue(residual_ji.jacobian_i[iVar][jVar]);
          Block_ji[iVar*nVar + jVar] -= 0.5 * SU2_TYPE::GetValue(residual_ji.jacobian_j[iVar][jVar]);
        }
    }
  } else {
    LinSysRes.SubtractBlock(jPoint, residual_ji);
    if (implicit) {
      /*--- The order of arguments were flipped in the evaluation of residual_ji, the Jacobian
       * associated with point i is stored in jacobian_j and point j in jacobian_i. ---*/
      for (int iVar=0; iVar<nVar; iVar++)
        for (int jVar=0; jVar<nVar; jVar++) {
          Block_ji[iVar*nVar + jVar] -= SU2_TYPE::GetValue(residual_ji.jacobian_j[iVar][jVar]);
          Block_jj[iVar*nVar + jVar] -= SU2_TYPE::GetValue(residual_ji.jacobian_i[iVar][jVar]);
        }
    }
  }
}
