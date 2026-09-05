/*!
 * \file nemo_viscous_assembly.cpp
 * \brief Unit test for the edge assembly of the NEMO viscous residual and Jacobian.
 * \author J. Bellon
 * \version 8.5.0 "Harrier"
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

#include "catch.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>
#include <vector>
#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../SU2_CFD/include/numerics/NEMO/NEMO_diffusion.hpp"
#include "../../SU2_CFD/include/solvers/CNEMONSSolver.hpp"
#include "../../SU2_CFD/include/variables/CNEMONSVariable.hpp"

namespace {

/*!
 * \brief A NEMO Navier-Stokes solver on a small box mesh, built the same way
 * the gradient tests build their geometry, with a non-uniform velocity field
 * so that the viscous residual is not trivially zero.
 */
struct NEMOViscousAssemblyCase {
  const std::string configOptions =
      "SOLVER= NEMO_NAVIER_STOKES\n"
      "GAS_MODEL= AIR-5\n"
      "GAS_COMPOSITION= (0.77, 0.23, 0.0, 0.0, 0.0)\n"
      "FLUID_MODEL= SU2_NONEQ\n"
      "MESH_FORMAT= BOX\n"
      "MESH_BOX_SIZE= 4,4,4\n"
      "MESH_BOX_LENGTH= 1,1,1\n"
      "MESH_BOX_OFFSET= 0,0,0\n"
      "INIT_OPTION= TD_CONDITIONS\n"
      "MACH_NUMBER= 0.3\n"
      "FREESTREAM_PRESSURE= 1000.0\n"
      "FREESTREAM_TEMPERATURE= 300.0\n"
      "FREESTREAM_TEMPERATURE_VE= 300.0\n"
      "REYNOLDS_NUMBER= 1000\n"
      "KIND_TURB_MODEL= NONE\n"
      "MARKER_FAR= (x_minus, x_plus, y_minus, y_plus, z_minus, z_plus)\n"
      "NUM_METHOD_GRAD= GREEN_GAUSS\n"
      "CONV_NUM_METHOD_FLOW= AUSM\n"
      "MUSCL_FLOW= NO\n"
      "TIME_DISCRE_FLOW= EULER_IMPLICIT\n";

  std::unique_ptr<CConfig> config;
  std::unique_ptr<CGeometry> geometry;
  CNEMONSSolver* solver{nullptr};

  NEMOViscousAssemblyCase() {
    auto origBuf = cout.rdbuf();
    cout.rdbuf(nullptr);

    stringstream ss(configOptions);
    config = std::unique_ptr<CConfig>(new CConfig(ss, SU2_COMPONENT::SU2_CFD, false));

    {
      auto aux_geometry = std::unique_ptr<CGeometry>(new CPhysicalGeometry(config.get(), 0, 1));
      geometry = std::unique_ptr<CGeometry>(new CPhysicalGeometry(aux_geometry.get(), config.get()));
    }
    geometry->SetSendReceive(config.get());
    geometry->SetBoundaries(config.get());
    geometry->SetPoint_Connectivity();
    geometry->SetElement_Connectivity();
    geometry->SetBoundVolume();
    geometry->Check_IntElem_Orientation(config.get());
    geometry->Check_BoundElem_Orientation(config.get());
    geometry->SetEdges();
    geometry->SetVertex(config.get());
    geometry->SetControlVolume(config.get(), ALLOCATE);
    geometry->SetBoundControlVolume(config.get(), ALLOCATE);
    geometry->FindNormal_Neighbor(config.get());
    geometry->SetGlobal_to_Local_Point();
    geometry->PreprocessP2PComms(geometry.get(), config.get());

    solver = new CNEMONSSolver(geometry.get(), config.get(), MESH_0);

    /*--- Scale the momentum with position so that the velocity gradients,
     * and with them the viscous fluxes, are non-zero. ---*/
    const auto nDim = geometry->GetnDim();
    const auto nSpecies = config->GetnSpecies();
    auto* nodes = solver->GetNodes();
    for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); ++iPoint) {
      const auto* coord = geometry->nodes->GetCoord(iPoint);
      const su2double scale = 1.0 + 0.2 * coord[0] + 0.1 * coord[1];
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
        const auto iVar = nSpecies + iDim;
        nodes->SetSolution(iPoint, iVar, scale * nodes->GetSolution(iPoint, iVar));
      }
    }

    cout.rdbuf(origBuf);
  }

  ~NEMOViscousAssemblyCase() { delete solver; }
};

}  // namespace

TEST_CASE("NEMO viscous edge assembly is the derivative of the subtracted and added residual",
          "[NEMO][viscous][Jacobian]") {
  NEMOViscousAssemblyCase testCase;
  auto* config = testCase.config.get();
  auto* geometry = testCase.geometry.get();
  auto* solver = testCase.solver;

  const auto nDim = geometry->GetnDim();
  const auto nVar = solver->GetnVar();
  const auto nPrimVar = solver->GetnPrimVar();
  const auto nPrimVarGrad = solver->GetnPrimVarGrad();
  const auto nPoint = geometry->GetnPoint();

  /*--- Primitive variables, transport properties and gradients. ---*/
  CSolver* solver_container[MAX_SOLS] = {nullptr};
  solver_container[FLOW_SOL] = solver;
  {
    auto origBuf = cout.rdbuf();
    cout.rdbuf(nullptr);
    solver->Preprocessing(geometry, solver_container, config, MESH_0, 0, RUNTIME_FLOW_SYS, false);
    cout.rdbuf(origBuf);
  }
  solver->LinSysRes.SetValZero();
  solver->Jacobian.SetValZero();

  /*--- Assemble the viscous residual and Jacobian with the solver. The
   * solver's override is private, so dispatch through the public base
   * interface, exactly as the integration classes do. ---*/
  CAvgGradCorrected_NEMO numerics(nDim, nVar, nPrimVar, nPrimVarGrad, config);
  CNumerics* numerics_container[MAX_TERMS] = {nullptr};
  numerics_container[VISC_TERM] = &numerics;
  CSolver& base = *solver;
  base.Viscous_Residual(geometry, solver_container, numerics_container, config, MESH_0, 0);

  /*--- Replay every edge with an independent numerics object. The residual
   * contribution F of an edge is subtracted at i and added at j, so the
   * derivative of the assembled residual is -dF/dU in the i rows and +dF/dU
   * in the j rows, with dF/dU_i and dF/dU_j the two blocks returned by the
   * numerics. Accumulate that reference and compare it with what the solver
   * assembled. ---*/
  CAvgGradCorrected_NEMO replay(nDim, nVar, nPrimVar, nPrimVarGrad, config);
  auto* nodes = dynamic_cast<CNEMONSVariable*>(solver->GetNodes());
  REQUIRE(nodes != nullptr);

  std::vector<su2double> residual_ref(nPoint * nVar, 0.0);
  std::vector<su2double> diagonal_ref(nPoint * nVar * nVar, 0.0);
  su2double offdiag_error = 0.0, jacobian_scale = 0.0;

  for (auto iEdge = 0ul; iEdge < geometry->GetnEdge(); ++iEdge) {
    const auto iPoint = geometry->edges->GetNode(iEdge, 0);
    const auto jPoint = geometry->edges->GetNode(iEdge, 1);

    replay.SetCoord(geometry->nodes->GetCoord(iPoint), geometry->nodes->GetCoord(jPoint));
    replay.SetNormal(geometry->edges->GetNormal(iEdge));
    replay.SetConservative(nodes->GetSolution(iPoint), nodes->GetSolution(jPoint));
    replay.SetPrimitive(nodes->GetPrimitive(iPoint), nodes->GetPrimitive(jPoint));
    replay.SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nodes->GetGradient_Primitive(jPoint));
    replay.SetdPdU(nodes->GetdPdU(iPoint), nodes->GetdPdU(jPoint));
    replay.SetdTdU(nodes->GetdTdU(iPoint), nodes->GetdTdU(jPoint));
    replay.SetdTvedU(nodes->GetdTvedU(iPoint), nodes->GetdTvedU(jPoint));
    replay.SetEve(nodes->GetEve(iPoint), nodes->GetEve(jPoint));
    replay.SetCvve(nodes->GetCvve(iPoint), nodes->GetCvve(jPoint));
    replay.SetDiffusionCoeff(nodes->GetDiffusionCoeff(iPoint), nodes->GetDiffusionCoeff(jPoint));
    replay.SetLaminarViscosity(nodes->GetLaminarViscosity(iPoint), nodes->GetLaminarViscosity(jPoint));
    replay.SetEddyViscosity(nodes->GetEddyViscosity(iPoint), nodes->GetEddyViscosity(jPoint));
    replay.SetThermalConductivity(nodes->GetThermalConductivity(iPoint), nodes->GetThermalConductivity(jPoint));
    replay.SetThermalConductivity_ve(nodes->GetThermalConductivity_ve(iPoint),
                                     nodes->GetThermalConductivity_ve(jPoint));

    const auto edge = replay.ComputeResidual(config);

    const auto* block_ij = solver->Jacobian.GetBlock(iPoint, jPoint);
    const auto* block_ji = solver->Jacobian.GetBlock(jPoint, iPoint);
    REQUIRE(block_ij != nullptr);
    REQUIRE(block_ji != nullptr);

    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      residual_ref[iPoint * nVar + iVar] -= edge.residual[iVar];
      residual_ref[jPoint * nVar + iVar] += edge.residual[iVar];

      for (unsigned short jVar = 0; jVar < nVar; ++jVar) {
        const su2double dFdUi = edge.jacobian_i[iVar][jVar];
        const su2double dFdUj = edge.jacobian_j[iVar][jVar];
        jacobian_scale = std::max(jacobian_scale, std::max(std::fabs(dFdUi), std::fabs(dFdUj)));

        /*--- Diagonal blocks collect every edge of a point. ---*/
        diagonal_ref[(iPoint * nVar + iVar) * nVar + jVar] -= dFdUi;
        diagonal_ref[(jPoint * nVar + iVar) * nVar + jVar] += dFdUj;

        /*--- Off-diagonal blocks belong to this edge alone. ---*/
        const su2double assembled_ij = block_ij[iVar * nVar + jVar];
        const su2double assembled_ji = block_ji[iVar * nVar + jVar];
        offdiag_error = std::max(offdiag_error, std::fabs(assembled_ij - (-dFdUj)));
        offdiag_error = std::max(offdiag_error, std::fabs(assembled_ji - (+dFdUi)));
      }
    }
  }

  su2double diagonal_error = 0.0, residual_error = 0.0, residual_scale = 0.0;
  for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
    const auto* block_ii = solver->Jacobian.GetBlock(iPoint, iPoint);
    const auto* assembled_residual = solver->LinSysRes.GetBlock(iPoint);
    REQUIRE(block_ii != nullptr);
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      residual_scale = std::max(residual_scale, std::fabs(residual_ref[iPoint * nVar + iVar]));
      residual_error =
          std::max(residual_error, std::fabs(assembled_residual[iVar] - residual_ref[iPoint * nVar + iVar]));
      for (unsigned short jVar = 0; jVar < nVar; ++jVar) {
        const su2double assembled = block_ii[iVar * nVar + jVar];
        diagonal_error =
            std::max(diagonal_error, std::fabs(assembled - diagonal_ref[(iPoint * nVar + iVar) * nVar + jVar]));
      }
    }
  }

  /*--- The reference must be non-trivial, otherwise the sign is untested. ---*/
  REQUIRE(jacobian_scale > 0.0);
  REQUIRE(residual_scale > 0.0);

  const su2double tolerance = 1.0e-10;
  CHECK(residual_error <= tolerance * residual_scale);
  CHECK(offdiag_error <= tolerance * jacobian_scale);
  CHECK(diagonal_error <= tolerance * jacobian_scale);
}
