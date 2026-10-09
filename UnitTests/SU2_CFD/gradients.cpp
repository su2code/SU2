/*!
 * \file gradients.cpp
 * \brief Unit tests for gradient calculation.
 * \author P. Gomes, T. Albring
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
#include "../../SU2_CFD/include/variables/CSpeciesFlameletVariable.hpp"
#include "../../SU2_CFD/include/solvers/CIncNSSolver.hpp"
#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../Common/include/containers/container_decorators.hpp"
#include "../../SU2_CFD/include/solvers/CSolver.hpp"
#include "../../SU2_CFD/include/gradients/computeGradientsGreenGauss.hpp"
#include "../../SU2_CFD/include/gradients/computeGradientsLeastSquares.hpp"
#include "../UnitQuadTestCase.hpp"

/*!
 * \brief Base class for gradient tests using a unit cube geometry.
 * Derived classes should implement operator (i,j), returning the value
 * of the test function, and method grad(i,j,k), returning the known
 * gradient.
 */
struct GradientTestBase {
  const string configOptions =
      "SOLVER= NAVIER_STOKES\n"
      "MESH_FORMAT= BOX\n"
      "INIT_OPTION= TD_CONDITIONS\n"
      "MARKER_HEATFLUX= (y_minus, 0.0, y_plus, 0.0)\n"
      "MARKER_FAR= (x_minus, x_plus, z_plus, z_minus)\n"
      "MESH_BOX_SIZE= 10,10,10\n"
      "MESH_BOX_LENGTH= 1,1,1\n"
      "MESH_BOX_OFFSET= 0,0,0\n";

  std::unique_ptr<CConfig> config;
  std::unique_ptr<CGeometry> geometry;

  GradientTestBase() {
    initConfig();
    initGeometry();
  }

  /*!
   * \brief Initialize the config structure
   */
  void initConfig() {
    auto origBuf = cout.rdbuf();
    cout.rdbuf(nullptr);
    stringstream ss(configOptions);
    config = std::unique_ptr<CConfig>(new CConfig(ss, SU2_COMPONENT::SU2_CFD, false));
    cout.rdbuf(origBuf);
  }

  /*!
   * \brief Initialize the geometry
   */
  void initGeometry() {
    auto origBuf = cout.rdbuf();
    cout.rdbuf(nullptr);
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

    cout.rdbuf(origBuf);
  }
};

struct LinearFunction : public GradientTestBase {
  const unsigned long nVar = 1;
  const su2double constant = -1.0;
  const su2double slope[3] = {1.0, 2.0, 3.0};

  /*!
   * \brief Return manufactured value.
   */
  su2double operator()(unsigned long iPoint, unsigned long) const {
    const auto coord = geometry->nodes->GetCoord(iPoint);
    return constant + GeometryToolbox::DotProduct(geometry->GetnDim(), slope, coord);
  }

  /*!
   * \brief Return reference value.
   */
  su2double grad(unsigned long, unsigned long, unsigned long iDim) const { return slope[iDim]; }
};

template <class T, class U>
void check(const T& ref, const U& calc, su2double tol = 1e-9) {
  su2double err = 0.0;
  for (auto iPoint = 0ul; iPoint < calc.length(); ++iPoint) {
    for (auto iVar = 0ul; iVar < calc.rows(); ++iVar)
      for (auto iDim = 0ul; iDim < calc.cols(); ++iDim)
        err = max(err, abs(calc(iPoint, iVar, iDim) - ref.grad(iPoint, iVar, iDim)));
  }
  CHECK(err < tol);
}

template <class TestField>
void testGreenGauss() {
  TestField field;
  C3DDoubleMatrix gradient(field.geometry->GetnPoint(), field.nVar, field.geometry->GetnDim());

  computeGradientsGreenGauss(nullptr, MPI_QUANTITIES::SOLUTION, PERIODIC_NONE, *field.geometry.get(),
                             *field.config.get(), field, 0, field.nVar, -1, gradient);
  check(field, gradient);
}

template <class TestField>
void testLeastSquares(bool weighted) {
  TestField field;
  const auto nDim = field.geometry->GetnDim();
  C3DDoubleMatrix R(field.geometry->GetnPoint(), nDim, nDim);
  C3DDoubleMatrix gradient(field.geometry->GetnPoint(), field.nVar, nDim);

  computeGradientsLeastSquares(nullptr, MPI_QUANTITIES::SOLUTION, PERIODIC_NONE, *field.geometry.get(),
                               *field.config.get(), weighted, field, 0, field.nVar, -1, gradient, R);
  check(field, gradient);
}

TEST_CASE("GG", "[Gradients]") { testGreenGauss<LinearFunction>(); }

TEST_CASE("LS", "[Gradients]") { testLeastSquares<LinearFunction>(false); }

TEST_CASE("WLS", "[Gradients]") { testLeastSquares<LinearFunction>(true); }

TEST_CASE("Periodic auxiliary gradients", "[Gradients][Periodic]") {
  const bool rotation = GENERATE(false, true);
  const bool stretched = GENERATE(false, true);
  const auto method = GENERATE(std::string("GREEN_GAUSS"), std::string("WEIGHTED_LEAST_SQUARES"));
  UnitQuadTestCase field;
  const string custom = "MARKER_CUSTOM= ( x_minus, x_plus, z_plus, z_minus)";
  field.config_options.replace(field.config_options.find(custom), custom.size(), "MARKER_CUSTOM= (z_plus, z_minus)");
  field.AddOption(rotation ? "MARKER_PERIODIC= (x_minus,x_plus, 0,0.5,0.5, 90,0,0, 1,0,0)"
                           : "MARKER_PERIODIC= (x_minus,x_plus, 0,0,0, 0,0,0, 1,0,0)");
  field.AddOption("VORTICITY_CONFINEMENT= YES");
  field.AddOption("NUM_METHOD_GRAD= " + method);
  field.InitConfig();
  field.InitGeometry(true);
  if (stretched) {
    for (auto i = 0ul; i < field.geometry->GetnPoint(); ++i) {
      const auto x = field.geometry->nodes->GetCoord(i, 0);
      field.geometry->nodes->SetCoord(i, 0, x * x);
    }
    field.geometry->SetControlVolume(field.config.get(), UPDATE);
    field.geometry->SetBoundControlVolume(field.config.get(), UPDATE);
  }
  field.geometry->MatchPeriodic(field.config.get(), 1);
  field.geometry->PreprocessPeriodicComms(field.geometry.get(), field.config.get());
  field.InitSolver();

  auto* solver = field.solver[FLOW_SOL];
  auto* nodes = solver->GetNodes();
  REQUIRE(nodes->GetnAuxVar() == 1);
  for (auto iPoint = 0ul; iPoint < field.geometry->GetnPoint(); ++iPoint) {
    const auto* coord = field.geometry->nodes->GetCoord(iPoint);
    const auto y = coord[1] - 0.5;
    const auto z = coord[2] - 0.5;
    nodes->SetAuxVar(iPoint, 0,
                     (rotation ? 2.0 + y * y + z * z : 2.0 + 3.0 * coord[1] - coord[2]) + coord[0] * (1 - coord[0]));
  }
  if (method == "GREEN_GAUSS") {
    SU2_OMP_PARALLEL { solver->SetAuxVar_Gradient_GG(field.geometry.get(), field.config.get()); }
    END_SU2_OMP_PARALLEL
  } else {
    SU2_OMP_PARALLEL { solver->SetAuxVar_Gradient_LS(field.geometry.get(), field.config.get()); }
    END_SU2_OMP_PARALLEL
  }

  su2double error = 0.0;
  for (auto iPoint = 0ul; iPoint < field.geometry->GetnPointDomain(); ++iPoint) {
    const auto* coord = field.geometry->nodes->GetCoord(iPoint);
    if (rotation && (coord[1] == 0.0 || coord[1] == 1.0 || coord[2] == 0.0 || coord[2] == 1.0)) continue;
    const auto index = static_cast<unsigned short>(
        std::round((stretched ? sqrt(SU2_TYPE::GetValue(coord[0])) : SU2_TYPE::GetValue(coord[0])) * 4));
    auto xcoord = [&](unsigned short k) {
      const auto x = 0.25 * k;
      return stretched ? x * x : x;
    };
    const auto left = xcoord(index == 0 ? 3 : index - 1), right = xcoord(index == 4 ? 1 : index + 1);
    const auto hm = coord[0] - (index == 0 ? left - 1 : left);
    const auto hp = (index == 4 ? right + 1 : right) - coord[0];
    const auto fm = left * (1 - left), fp = right * (1 - right), fi = coord[0] * (1 - coord[0]);
    const auto gx = method == "GREEN_GAUSS" ? (fp - fm) / (hp + hm) : ((fp - fi) / hp + (fi - fm) / hm) / 2;
    const su2double expected[3] = {gx, rotation ? 2.0 * (coord[1] - 0.5) : 3.0,
                                   rotation ? 2.0 * (coord[2] - 0.5) : -1.0};
    for (auto iDim = 0u; iDim < 3; ++iDim)
      error = max(error, abs(nodes->GetAuxVarGradient(iPoint, 0, iDim) - expected[iDim]));
  }
  CHECK(error < 1e-9);
}

TEST_CASE("Flamelet preferential diffusion auxiliary gradients", "[Gradients][Flamelet][Periodic]") {
  const bool periodic = GENERATE(false, true);
  const auto method = GENERATE(std::string("GREEN_GAUSS"), std::string("WEIGHTED_LEAST_SQUARES"));
  UnitQuadTestCase field;
  field.SetOption("SOLVER= INC_NAVIER_STOKES");
  field.SetOption("KIND_VERIFICATION_SOLUTION= NO_VERIFICATION_SOLUTION");
  field.AddOption("KIND_SCALAR_MODEL= FLAMELET");
  field.AddOption("FLUID_MODEL= FLUID_FLAMELET");
  field.AddOption("INC_DENSITY_MODEL= VARIABLE");
  field.SetOption("VISCOSITY_MODEL= FLAMELET");
  field.AddOption("CONDUCTIVITY_MODEL= FLAMELET");
  field.AddOption("DIFFUSIVITY_MODEL= FLAMELET");
  field.AddOption(
      "SPECIES_INIT= (0,300000,0.5)\nCONTROLLING_VARIABLE_NAMES= (ProgressVariable,EnthalpyTot,MixtureFraction)");
  field.AddOption("CONTROLLING_VARIABLE_SOURCE_NAMES= (ProdRateTot_PV,NULL,NULL)");
  field.AddOption("PREFERENTIAL_DIFFUSION= YES\nNUM_METHOD_GRAD= " + method);
  if (periodic) {
    const string custom = "MARKER_CUSTOM= ( x_minus, x_plus, z_plus, z_minus)";
    field.config_options.replace(field.config_options.find(custom), custom.size(), "MARKER_CUSTOM= (z_plus,z_minus)");
    field.AddOption("MARKER_PERIODIC= (x_minus,x_plus, 0,0,0, 0,0,0, 1,0,0)");
  }
  field.InitConfig();
  field.InitGeometry(true);
  if (periodic) {
    field.geometry->MatchPeriodic(field.config.get(), 1);
    field.geometry->PreprocessPeriodicComms(field.geometry.get(), field.config.get());
  }
  const su2double initial[] = {0, 300000, 0.5};
  CSpeciesFlameletVariable nodes(initial, field.geometry->GetnPoint(), 3, 3, field.config.get());
  CHECK(nodes.GetnAuxVar() == FLAMELET_PREF_DIFF_SCALARS::N_BETA_TERMS);
  /*--- Exercise the shared auxiliary communication with four beta fields and
   * three transported scalars, without requiring an unrelated flamelet table. ---*/
  struct AuxiliarySolver : CSolver {
    CVariable* nodes;
    CVariable* GetBaseClassPointerToNodes() override { return nodes; }
    AuxiliarySolver(CVariable& values, CGeometry& geometry) : nodes(&values) {
      nDim = 3;
      nVar = 3;
      nPoint = geometry.GetnPoint();
      nPointDomain = geometry.GetnPointDomain();
      SetBaseClassPointerToNodes();
    }
  } solver(nodes, *field.geometry);
  if (periodic) {
    /*--- The flow solver normally supplies the complete periodic volume. ---*/
    SU2_OMP_PARALLEL {
      solver.InitiatePeriodicComms(field.geometry.get(), field.config.get(), 1, PERIODIC_VOLUME);
      solver.CompletePeriodicComms(field.geometry.get(), field.config.get(), 1, PERIODIC_VOLUME);
    }
    END_SU2_OMP_PARALLEL
  }
  for (auto i = 0ul; i < field.geometry->GetnPoint(); ++i) {
    const auto* x = field.geometry->nodes->GetCoord(i);
    for (auto v = 0u; v < FLAMELET_PREF_DIFF_SCALARS::N_BETA_TERMS; ++v)
      nodes.SetAuxVar(i, v, (v + 1) * (2 + (periodic ? x[0] * (1 - x[0]) : x[0]) + 3 * x[1] - x[2]));
  }
  if (method == "GREEN_GAUSS") {
    SU2_OMP_PARALLEL { solver.SetAuxVar_Gradient_GG(field.geometry.get(), field.config.get()); }
    END_SU2_OMP_PARALLEL
  } else {
    SU2_OMP_PARALLEL { solver.SetAuxVar_Gradient_LS(field.geometry.get(), field.config.get()); }
    END_SU2_OMP_PARALLEL
  }
  su2double error = 0;
  for (auto i = 0ul; i < field.geometry->GetnPointDomain(); ++i) {
    const auto* x = field.geometry->nodes->GetCoord(i);
    const su2double expected[] = {periodic ? (x[0] == 0 || x[0] == 1 ? 0.0 : 1 - 2 * x[0]) : 1.0, 3, -1};
    for (auto v = 0u; v < FLAMELET_PREF_DIFF_SCALARS::N_BETA_TERMS; ++v)
      for (auto d = 0u; d < 3; ++d)
        error = std::max(error, fabs(nodes.GetAuxVarGradient(i, v, d) - (v + 1) * expected[d]));
  }
  CHECK(error < 1e-12);
}

TEST_CASE("Recovered periodic wall heat flux", "[Periodic][Heat]") {
  UnitQuadTestCase field;
  field.SetOption("SOLVER= INC_NAVIER_STOKES");
  field.AddOption("FLUID_MODEL= CONSTANT_DENSITY");
  field.SetOption("KIND_VERIFICATION_SOLUTION= NO_VERIFICATION_SOLUTION");
  field.AddOption("INC_NONDIM= DIMENSIONAL");
  field.AddOption("INC_ENERGY_EQUATION= YES\nKIND_STREAMWISE_PERIODIC= PRESSURE_DROP");
  field.AddOption("STREAMWISE_PERIODIC_TEMPERATURE= YES\nINC_VELOCITY_INIT= (1,0,0)");
  const string walls = "MARKER_HEATFLUX= (y_minus, 0.0, y_plus, 0.0)";
  field.config_options.replace(field.config_options.find(walls), walls.size(), "MARKER_HEATFLUX= (y_minus,1,y_plus,3)");
  const string custom = "MARKER_CUSTOM= ( x_minus, x_plus, z_plus, z_minus)";
  field.config_options.replace(field.config_options.find(custom), custom.size(), "MARKER_CUSTOM= (z_plus,z_minus)");
  field.AddOption("MARKER_PERIODIC= (x_minus,x_plus, 0,0,0, 0,0,0, 1,0,0)");
  field.InitConfig();
  field.InitGeometry(true);
  for (auto i = 0ul; i < field.geometry->GetnPoint(); ++i) {
    const auto* x = field.geometry->nodes->GetCoord(i);
    field.geometry->nodes->SetCoord(i, 1, x[1] * (1 + 0.1 * sin(2 * PI_NUMBER * x[0])));
  }
  field.geometry->SetControlVolume(field.config.get(), UPDATE);
  field.geometry->SetBoundControlVolume(field.config.get(), UPDATE);
  field.geometry->MatchPeriodic(field.config.get(), 1);
  field.geometry->PreprocessPeriodicComms(field.geometry.get(), field.config.get());
  field.InitSolver();
  auto* solver = static_cast<CIncNSSolver*>(field.solver[FLOW_SOL]);
  field.config->SetGlobalParam(field.config->GetKind_Solver(), RUNTIME_FLOW_SYS);
  {
    SU2_OMP_PARALLEL {
      solver->Preprocessing(field.geometry.get(), field.solver, field.config.get(), MESH_0, 0, RUNTIME_FLOW_SYS, false);
    }
    END_SU2_OMP_PARALLEL
  }
  const auto values = solver->GetStreamwisePeriodicValues();
  REQUIRE(fabs(values.Streamwise_Periodic_MassFlow) > 1e-12);
  REQUIRE(fabs(values.Streamwise_Periodic_IntegratedHeatFlow) > 1e-12);
  auto* nodes = solver->GetNodes();
  su2double error = 0;
  unsigned long checked = 0, global = 0;
  for (auto m = 0u; m < field.geometry->GetnMarker(); ++m) {
    if (field.config->GetMarker_All_TagBound(m) != "y_plus") continue;
    solver->LinSysRes = su2double(0);
    {
      SU2_OMP_PARALLEL {
        solver->BC_HeatFlux_Wall(field.geometry.get(), field.solver, nullptr, nullptr, field.config.get(), m);
      }
      END_SU2_OMP_PARALLEL
    }
    for (auto v = 0ul; v < field.geometry->GetnVertex(m); ++v) {
      const auto i = field.geometry->vertex[m][v]->GetNode();
      if (!field.geometry->nodes->GetDomain(i)) continue;
      const auto* normal = field.geometry->vertex[m][v]->GetNormal();
      if (fabs(normal[0]) < 1e-12) continue;
      const auto area = GeometryToolbox::Norm(3, normal);
      /*--- T = T_periodic + g.x. With an inward area normal, the periodic
       * conductive heat flux entering the fluid is q + k g.n. ---*/
      const auto gx = values.Streamwise_Periodic_IntegratedHeatFlow /
                      (values.Streamwise_Periodic_MassFlow * nodes->GetSpecificHeatCp(i));
      const auto expected = -field.config->GetWall_HeatFlux("y_plus") / field.config->GetHeat_Flux_Ref() * area -
                            nodes->GetThermalConductivity(i) * gx * normal[0];
      error = std::max(error, fabs(solver->LinSysRes(i, 4) - expected));
      ++checked;
    }
  }
  SU2_MPI::Allreduce(&checked, &global, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  REQUIRE(global > 0);
  CHECK(error < 1e-12);
}
