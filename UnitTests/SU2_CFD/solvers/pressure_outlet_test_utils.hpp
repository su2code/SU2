/*!
 * \file pressure_outlet_test_utils.hpp
 * \brief Shared fixture and independent references for the static-pressure outlet tests.
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

#pragma once

#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include "catch.hpp"
#include "../../UnitQuadTestCase.hpp"
#include "../../../SU2_CFD/include/solvers/CEulerSolver.hpp"
#include "../../../SU2_CFD/include/variables/CNSVariable.hpp"
#include "../../../SU2_CFD/include/fluid/CIdealGas.hpp"
#include "../../../SU2_CFD/include/fluid/CPengRobinson.hpp"
#include "../../../SU2_CFD/include/fluid/CVanDerWaalsGas.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

/*--- Reference state of the tests: siloxane MDM near its critical point, or an ideal gas at the same p and rho. ---*/
constexpr double OUTLET_TEST_P = 1500112, OUTLET_TEST_RHO = 230;
constexpr unsigned short OUTLET_TEST_MAXVAR = 5;

// Gases of the outlet tests: 0 ideal (gamma 1.4), 1 Peng-Robinson, 2 van der Waals.
inline std::unique_ptr<CFluidModel> MakeOutletTestGas(unsigned short kind) {
  if (kind == 0) return std::unique_ptr<CFluidModel>(new CIdealGas(1.4, 35.152));
  if (kind == 1) return std::unique_ptr<CFluidModel>(new CPengRobinson(1.0125, 35.152, 1415200, 564.1, .529));
  return std::unique_ptr<CFluidModel>(new CVanDerWaalsGas(1.0125, 35.152, 1415200, 564.1));
}

// Configuration of a box case (x_minus..z_plus) with the gas of MakeOutletTestGas; the caller adds the markers.
inline std::string OutletTestOptions(unsigned short dim, unsigned short kind, double turbulentK = 0,
                                     bool implicit = true, const std::string& boxSize = "5,5,5") {
  std::ostringstream options;
  options << std::setprecision(17) << "SOLVER= " << (turbulentK ? "RANS" : "EULER")
          << "\nMESH_FORMAT= " << (dim == 2 ? "RECTANGLE" : "BOX") << "\nMESH_BOX_SIZE= " << boxSize
          << "\nREF_DIMENSIONALIZATION= DIMENSIONAL\nMACH_NUMBER= 0.3\nFREESTREAM_OPTION= DENSITY_FS"
          << "\nFREESTREAM_DENSITY= " << OUTLET_TEST_RHO << "\nFREESTREAM_PRESSURE= " << OUTLET_TEST_P
          << "\nFLUID_MODEL= "
          << (kind == 0   ? "IDEAL_GAS"
              : kind == 1 ? "PR_GAS"
                          : "VW_GAS")
          << "\nGAMMA_VALUE= " << (kind == 0 ? "1.4" : "1.0125")
          << "\nGAS_CONSTANT= 35.152\nCRITICAL_PRESSURE= 1415200\nCRITICAL_TEMPERATURE= 564.1\nACENTRIC_FACTOR= 0.529"
          << "\nCONV_NUM_METHOD_FLOW= HLLC\nTIME_DISCRE_FLOW= "
          << (implicit ? "EULER_IMPLICIT" : "RUNGE-KUTTA_EXPLICIT") << '\n';
  if (turbulentK)
    options << "KIND_TURB_MODEL= SST\nINIT_OPTION= TD_CONDITIONS\nREYNOLDS_NUMBER= 1E6\n"
               "VISCOSITY_MODEL= CONSTANT_VISCOSITY\nCONDUCTIVITY_MODEL= CONSTANT_CONDUCTIVITY\n";
  return options.str();
}

// Conservative state with the reference p and rho, the given velocity and turbulent kinetic energy.
inline void OutletTestState(CFluidModel& gas, unsigned short dim, const su2double* velocity, double turbulentK,
                            su2double* state) {
  gas.SetTDState_Prho(OUTLET_TEST_P, OUTLET_TEST_RHO);
  state[0] = OUTLET_TEST_RHO;
  for (auto d = 0u; d < dim; ++d) state[d + 1] = OUTLET_TEST_RHO * velocity[d];
  state[dim + 1] =
      OUTLET_TEST_RHO * (gas.GetStaticEnergy() + .5 * GeometryToolbox::SquaredNorm(dim, velocity) + turbulentK);
}

// Uniform flow (and grid velocity) at every point, with zero residual and Jacobian.
inline void SetOutletTestState(UnitQuadTestCase& test, CFluidModel& gas, const su2double* state, const su2double* grid,
                               double turbulentK) {
  if (turbulentK) {
    gas.SetLaminarViscosityModel(test.config.get());
    gas.SetThermalConductivityModel(test.config.get());
  }
  auto* solver = test.solver[FLOW_SOL];
  auto* nodes = static_cast<CEulerVariable*>(solver->GetNodes());
  const auto nVar = solver->GetnVar();
  for (auto point = 0ul; point < test.geometry->GetnPoint(); ++point) {
    nodes->SetSolution(point, state);
    nodes->SetSolution_Old(point, state);
    if (turbulentK) {
      test.solver[TURB_SOL]->GetNodes()->SetSolution(point, 0, turbulentK);
      REQUIRE(static_cast<CNSVariable*>(nodes)->SetPrimVar(point, 0, turbulentK, &gas));
    } else {
      REQUIRE(nodes->SetPrimVar(point, &gas));
    }
    nodes->SetSecondaryVar(point, &gas);
    if (grid != nullptr) test.geometry->nodes->SetGridVel(point, grid);
    for (auto k = 0u; k < nVar; ++k) solver->LinSysRes(point, k) = 0;
    su2double zero[OUTLET_TEST_MAXVAR * OUTLET_TEST_MAXVAR] = {};
    solver->Jacobian.SetBlock(point, point, zero);
  }
}

inline unsigned short OutletTestMarker(const CConfig& config, const std::string& tag) {
  for (auto marker = 0u; marker < config.GetnMarker_All(); ++marker)
    if (config.GetMarker_All_TagBound(marker) == tag) return marker;
  FAIL("marker " << tag << " not found");
  return 0;
}

// Outward area normal and unit normal of a boundary vertex.
inline void OutwardNormal(const CGeometry& geometry, unsigned short marker, unsigned long vertex, unsigned short dim,
                          su2double* normal, su2double* unit) {
  geometry.vertex[marker][vertex]->GetNormal(normal);
  const su2double area = GeometryToolbox::Norm(dim, normal);
  for (auto d = 0u; d < dim; ++d) unit[d] = (normal[d] = -normal[d]) / area;
}

// Mock numerics that return the boundary state as the residual and unit flux Jacobians, or zero (no diffusion).
class OutletFluxProbe final : public CNumerics {
  su2double residual[OUTLET_TEST_MAXVAR] = {}, matrix[OUTLET_TEST_MAXVAR][OUTLET_TEST_MAXVAR] = {};
  const su2double* rows[OUTLET_TEST_MAXVAR];
  bool zero;

 public:
  OutletFluxProbe(unsigned short dim, const CConfig* config, bool zeroFlux = false)
      : CNumerics(dim, dim + 2, config), zero(zeroFlux) {
    for (auto k = 0u; k < nVar; ++k) {
      rows[k] = matrix[k];
      matrix[k][k] = zero ? 0 : 1;
    }
  }
  ResidualType<> ComputeResidual(const CConfig*) override {
    if (!zero) {
      const su2double rho = V_j[nDim + 2], p = V_j[nDim + 1];
      residual[0] = rho;
      for (auto d = 0u; d < nDim; ++d) residual[d + 1] = rho * V_j[d + 1];
      residual[nDim + 1] = rho * V_j[nDim + 3] - p;
    }
    return {residual, rows, rows};
  }
};

// Physical boundary flux F(U_b).N - (v_grid.N) U_b of the Riemann outlets, whose EOS energy contains the SST k.
inline void PhysicalOutletFlux(CFluidModel& gas, const CNumerics& numerics, unsigned short dim,
                               const su2double* boundary, const su2double* normal, const su2double* grid,
                               su2double* flux) {
  su2double velocity[3];
  for (auto d = 0u; d < dim; ++d) velocity[d] = boundary[d + 1] / boundary[0];
  const su2double energy = boundary[dim + 1] / boundary[0];
  gas.SetTDState_rhoe(boundary[0], energy - .5 * GeometryToolbox::SquaredNorm(dim, velocity));
  const su2double pressure = gas.GetPressure(), enthalpy = energy + pressure / boundary[0];
  numerics.GetInviscidProjFlux(&boundary[0], velocity, &pressure, &enthalpy, normal, flux);
  const su2double gridFlux = grid ? GeometryToolbox::DotProduct(dim, grid, normal) : 0.0;
  for (auto k = 0u; k < dim + 2u; ++k) flux[k] -= gridFlux * boundary[k];
}

// Primitive boundary state that the scalar solvers read: it must match the EOS state of the last flux evaluation.
inline void CheckOutletPrimitives(const su2double* primitive, const CFluidModel& gas, unsigned short dim) {
  CHECK(SU2_TYPE::GetValue(primitive[0]) == Approx(SU2_TYPE::GetValue(gas.GetTemperature())).epsilon(1e-8));
  CHECK(SU2_TYPE::GetValue(primitive[dim + 1]) == Approx(SU2_TYPE::GetValue(gas.GetPressure())).epsilon(1e-8));
  CHECK(SU2_TYPE::GetValue(primitive[dim + 2]) == Approx(SU2_TYPE::GetValue(gas.GetDensity())).epsilon(1e-8));
  CHECK(SU2_TYPE::GetValue(primitive[dim + 4]) == Approx(SU2_TYPE::GetValue(gas.GetSoundSpeed())).epsilon(1e-8));
}

// Independent primitive-variable wave decomposition of the outgoing-wave projection.
inline void OutgoingOutletIncrement(CFluidModel& gas, unsigned short dim, const su2double* state,
                                    const su2double* normal, su2double gridNormal, su2double turbulentK,
                                    unsigned short column, su2double* increment) {
  const auto rho = state[0];
  su2double velocity[3] = {}, dv[3] = {};
  const su2double drho = column == 0 ? 1.0 : 0.0;
  for (auto d = 0u; d < dim; ++d) {
    velocity[d] = state[d + 1] / rho;
    dv[d] = ((column == d + 1 ? 1.0 : 0.0) - velocity[d] * drho) / rho;
  }
  const auto energy = state[dim + 1] / rho;
  gas.SetTDState_rhoe(rho, energy - .5 * GeometryToolbox::SquaredNorm(dim, velocity) - turbulentK);
  const auto sound = gas.GetSoundSpeed(), dpde = gas.GetdPde_rho(), dpdrho = gas.GetdPdrho_e();
  const su2double de =
      ((column == dim + 1 ? 1.0 : 0.0) - energy * drho) / rho - GeometryToolbox::DotProduct(dim, velocity, dv);
  const su2double dp = dpdrho * drho + dpde * de;
  const su2double q = GeometryToolbox::DotProduct(dim, velocity, normal) - gridNormal;
  const su2double dq = GeometryToolbox::DotProduct(dim, dv, normal);
  su2double projectedRho = 0, projectedP = 0, projectedVelocity[3] = {};
  if (q >= 0) {
    projectedRho = drho - dp / (sound * sound);
    for (auto d = 0u; d < dim; ++d) projectedVelocity[d] = dv[d] - dq * normal[d];
  }
  for (const auto sign : {-1, 1}) {
    if (q + sign * sound < 0) continue;
    const su2double amplitude = .5 * (dp + sign * rho * sound * dq);
    projectedRho += amplitude / (sound * sound);
    projectedP += amplitude;
    for (auto d = 0u; d < dim; ++d) projectedVelocity[d] += sign * amplitude / (rho * sound) * normal[d];
  }
  increment[0] = projectedRho;
  for (auto d = 0u; d < dim; ++d) increment[d + 1] = velocity[d] * projectedRho + rho * projectedVelocity[d];
  increment[dim + 1] = energy * projectedRho + rho * (projectedP - dpdrho * projectedRho) / dpde +
                       rho * GeometryToolbox::DotProduct(dim, velocity, projectedVelocity);
}

// Boundary state that the static-pressure outlets of a steady run use (CEulerSolver::BuildPressureOutletState):
// the outlet state, the sonic backflow state for flow entering through every characteristic, or else the interior.
// Returns true if a fallback was used.
inline bool ExpectedSteadyOutletState(CFluidModel& gas, unsigned short dim, const su2double* state,
                                      const su2double* normal, su2double gridNormal, su2double pressure,
                                      su2double turbulentK, bool idealGas, su2double gamma, su2double* boundary) {
  using Status = CEulerSolver::PressureOutletStatus;
  auto status = CEulerSolver::ComputePressureOutletState(&gas, dim, state, normal, gridNormal, pressure, turbulentK,
                                                         idealGas, gamma, boundary);
  if (status == Status::Valid) return false;
  if (status == Status::Backflow)
    status = CEulerSolver::ComputeSonicBackflowState(&gas, dim, state, normal, gridNormal, pressure, turbulentK,
                                                     idealGas, gamma, boundary);
  if (status != Status::Valid)
    for (auto k = 0u; k < dim + 2u; ++k) boundary[k] = state[k];
  return true;
}

// Fallback count of the last boundary call: the monitor sums and resets it (its warnings are discarded here).
inline unsigned long CountOutletFallbacks(CSolver* solver, const CConfig* config) {
  std::ostringstream quiet;
  auto* buffer = std::cout.rdbuf(quiet.rdbuf());
  auto* euler = static_cast<CEulerSolver*>(solver);
  euler->MonitorPressureOutletFallbacks(config);
  std::cout.rdbuf(buffer);
  return euler->GetPressureOutletFallbacks();
}

// Implicit block of a Riemann outlet: the physical flux Jacobian at the boundary state along the outgoing waves of the
// interior state, checked by finite differences. The block may be stored in reduced precision.
template <class Scalar>
void CheckRiemannOutletJacobian(const Scalar* block, CFluidModel& gas, const CNumerics& numerics, unsigned short dim,
                                const su2double* state, const su2double* boundary, const su2double* unit,
                                const su2double* normal, const su2double* grid) {
  const su2double gridNormal = grid ? GeometryToolbox::DotProduct(dim, grid, unit) : 0.0;
  gas.SetTDState_rhoe(
      state[0], state[dim + 1] / state[0] - .5 * GeometryToolbox::SquaredNorm(dim, state + 1) / (state[0] * state[0]));
  const double rho = SU2_TYPE::GetValue(state[0]), sound = SU2_TYPE::GetValue(gas.GetSoundSpeed());
  for (auto column = 0u; column < dim + 2u; ++column) {
    su2double plus[5], minus[5], fp[5], fm[5], increment[5];
    OutgoingOutletIncrement(gas, dim, state, unit, gridNormal, 0, column, increment);
    const double delta = 2e-5 * fmax(fabs(SU2_TYPE::GetValue(state[column])), column == 0 ? rho : rho * sound);
    for (auto k = 0u; k < dim + 2u; ++k) {
      plus[k] = boundary[k] + delta * increment[k];
      minus[k] = boundary[k] - delta * increment[k];
    }
    PhysicalOutletFlux(gas, numerics, dim, plus, normal, grid, fp);
    PhysicalOutletFlux(gas, numerics, dim, minus, normal, grid, fm);
    for (auto row = 0u; row < dim + 2u; ++row)
      CHECK(block[row * (dim + 2) + column] ==
            Approx(SU2_TYPE::GetValue((fp[row] - fm[row]) / (2 * delta))).epsilon(5e-4).margin(2e-5));
  }
}
