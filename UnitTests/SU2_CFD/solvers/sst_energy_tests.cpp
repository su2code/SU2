/*!
 * \file sst_energy_tests.cpp
 * \brief SST thermal-energy relaxation and matrix dissipation regressions.
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
#include "../../UnitQuadTestCase.hpp"

TEST_CASE("SST relaxation limits thermal energy when k dominates", "[SST][UnderRelaxation]") {
  UnitQuadTestCase test;
  test.config_options =
      "SOLVER= RANS\nMACH_NUMBER= 0.5\nREYNOLDS_NUMBER= 1e6\nKIND_TURB_MODEL= SST\nMESH_FORMAT= RECTANGLE\nMESH_BOX_SIZE= 5,5,5\n"
      "MARKER_HEATFLUX= (x_minus,0,x_plus,0,y_minus,0,y_plus,0)\n"
      "TIME_DISCRE_FLOW= EULER_IMPLICIT\nMAX_UPDATE_FLOW= 0.2\n";
  std::stringstream options(test.config_options);
  test.config = std::make_unique<CConfig>(options, SU2_COMPONENT::SU2_CFD, false);
  test.InitGeometry();
  test.solver = CSolverFactory::CreateSolverContainer(test.config->GetKind_Solver(), test.config.get(),
                                                      test.geometry.get(), 0);
  auto* flow = static_cast<CEulerSolver*>(test.solver[FLOW_SOL]);
  auto* nodes = flow->GetNodes();
  auto* turb = test.solver[TURB_SOL]->GetNodes();
  for (auto p = 0ul; p < test.geometry->GetnPoint(); ++p) {
    nodes->SetSolution(p, 0, 2.0);
    nodes->SetSolution(p, 1, 0.0);
    nodes->SetSolution(p, 2, 0.0);
    nodes->SetSolution(p, 3, 201.0);
    turb->SetSolution(p, 0, 100.0);
    for (auto v = 0u; v < 4; ++v) flow->LinSysSol(p, v) = 0.0;
    flow->LinSysSol(p, 3) = -10.0;
  }
  flow->ComputeUnderRelaxationFactor(test.config.get(), test.solver);
  CHECK(nodes->GetUnderRelaxation(0) == Approx(0.02));
  CHECK(201.0 - 10.0 * nodes->GetUnderRelaxation(0) - 200.0 > 0.0);
  // Exercise the production implicit-Euler caller, not only the relaxation helper.
  nodes->Set_OldSolution();
  flow->Jacobian.SetValZero();
  for (auto p = 0ul; p < test.geometry->GetnPoint(); ++p) {
    nodes->SetDelta_Time(p, test.geometry->nodes->GetVolume(p));
    nodes->SetRes_TruncErrorZero(p);
    for (auto v = 0u; v < 4; ++v) flow->LinSysRes(p, v) = 0.0;
    flow->LinSysRes(p, 3) = 10.0;
  }
  flow->ImplicitEuler_Iteration(test.geometry.get(), test.solver, test.config.get());
  CHECK(nodes->GetUnderRelaxation(0) == Approx(0.02));
  CHECK(nodes->GetSolution(0, 3) == Approx(200.8));
  delete test.solver[TURB_SOL];
}

#include "../../../SU2_CFD/include/variables/CNSVariable.hpp"
#include "../../../SU2_CFD/include/numerics_simd/flow/diffusion/viscous_fluxes.hpp"
#include "../../../SU2_CFD/include/numerics_simd/flow/convection/centered.hpp"

TEST_CASE("JST matrix fourth-order k dissipation does not create pressure waves", "[SST][JST]") {
  UnitQuadTestCase test;
  test.config_options =
      "SOLVER= RANS\nMACH_NUMBER= 0.5\nREYNOLDS_NUMBER= 1e6\nKIND_TURB_MODEL= SST\nMESH_FORMAT= RECTANGLE\nMESH_BOX_SIZE= 5,5,5\n"
      "MARKER_HEATFLUX= (x_minus,0,x_plus,0,y_minus,0,y_plus,0)\n"
      "CONV_NUM_METHOD_FLOW= JST_MAT\nENTROPY_FIX_COEFF= 0\n";
  std::stringstream options(test.config_options);
  test.config = std::make_unique<CConfig>(options, SU2_COMPONENT::SU2_CFD, false);
  test.InitGeometry();
  test.solver = CSolverFactory::CreateSolverContainer(test.config->GetKind_Solver(), test.config.get(),
                                                      test.geometry.get(), 0);
  auto* flow = static_cast<CEulerSolver*>(test.solver[FLOW_SOL]);
  auto& nodes = *static_cast<CEulerVariable*>(flow->GetNodes());
  const CVariable* turb = test.solver[TURB_SOL]->GetNodes();
  for (auto p = 0ul; p < test.geometry->GetnPoint(); ++p) {
    nodes.SetSensor(p, 0.0);
    for (auto v = 0u; v < 5; ++v) nodes.SetUnd_Lapl(p, v, 0.0);
  }
  nodes.SetUnd_Lapl(0, 3, 2.0);
  nodes.SetUnd_Lapl(0, 4, 2.0);  // the energy Laplacian comes entirely from rho*k
  CJSTmatScheme<CNoViscousFlux<2>> scheme(*test.config, turb);
  CPair<CCompressiblePrimitives<2, 7>> V;
  V.i.all = Double(0.0);
  V.i.density() = 1.0;
  V.i.pressure() = 1.0;
  V.i.speedSound() = sqrt(1.4);
  V.i.enthalpy() = 3.5;
  V.j = V.i;
  VectorDbl<4> flux, diff;
  flux = Double(0.0);
  diff = Double(0.0);
  VectorDbl<2> normal;
  normal(0) = 1.0;
  normal(1) = 0.0;
  MatrixDbl<4> ji, jj;
  scheme.finalizeFlux(flux, ji, jj, false, Double(1.0), Double(0.0), V.i, V, diff, Int(0ul), Int(1ul), *test.geometry,
                      nodes, normal);
  for (size_t lane = 0; lane < Double::Size; ++lane)
    for (size_t v = 0; v < 4; ++v) CHECK(flux(v)[lane] == Approx(0.0).margin(1e-12));
  delete test.solver[TURB_SOL];
}
