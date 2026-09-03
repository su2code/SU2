/*!
 * \file CNumerics_tests.cpp
 * \brief Unit tests for the numerics classes.
 * \author C. Pederson
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
#include <array>
#include <sstream>
#include <vector>
#include "../../../SU2_CFD/include/numerics/CNumerics.hpp"
#include "../../../SU2_CFD/include/numerics/NEMO/NEMO_diffusion.hpp"

TEST_CASE("NTS blending has a minimum of 0.05", "[Upwind/central blending]") {
  std::stringstream config_options;

  config_options << "SOLVER= NAVIER_STOKES" << std::endl;
  config_options << "ROE_LOW_DISSIPATION= "
                 << "NTS" << std::endl;
  config_options << "REYNOLDS_NUMBER= 5" << std::endl;

  /*--- Setup ---*/

  CConfig* config = new CConfig(config_options, SU2_COMPONENT::SU2_CFD, false);

  const su2double dissipation_i = 0;
  const su2double dissipation_j = 0;
  const su2double sensor_i = 0;
  const su2double sensor_j = 0;

  /*--- Test ---*/

  CNumerics numerics;
  const su2double dissipation_ij =
      numerics.GetRoe_Dissipation(dissipation_i, dissipation_j, sensor_i, sensor_j, config);

  REQUIRE(dissipation_ij >= 0.05);

  /*--- Teardown ---*/

  delete config;
}

TEST_CASE("QCR2000 corrects only the turbulent stress", "[QCR]") {
  constexpr size_t nDim = 3;
  su2double gradvel[3][3] = {{0.0}};
  gradvel[0][1] = 100.0;
  const su2double mu_lam = 2.0;
  const su2double mu_turb = 3.0;

  /*--- At a no-slip wall the eddy viscosity is zero, so the correction
   * must leave the (purely laminar) stress tensor unchanged. ---*/

  su2double tau_wall[3][3] = {{0.0}}, tau_wall_ref[3][3] = {{0.0}};
  CNumerics::ComputeStressTensor(nDim, tau_wall, gradvel, mu_lam);
  CNumerics::ComputeStressTensor(nDim, tau_wall_ref, gradvel, mu_lam);
  CNumerics::AddQCR(nDim, gradvel, tau_wall, 0.0);
  for (size_t iDim = 0; iDim < nDim; iDim++)
    for (size_t jDim = 0; jDim < nDim; jDim++) REQUIRE(tau_wall[iDim][jDim] == tau_wall_ref[iDim][jDim]);

  /*--- Scaling the correction of the total stress by mu_t / (mu_l + mu_t)
   * is equivalent to correcting the turbulent stress tensor alone. ---*/

  su2double tau_total[3][3] = {{0.0}}, tau_turb[3][3] = {{0.0}}, tau_lam[3][3] = {{0.0}};
  CNumerics::ComputeStressTensor(nDim, tau_total, gradvel, mu_lam + mu_turb);
  CNumerics::ComputeStressTensor(nDim, tau_turb, gradvel, mu_turb);
  CNumerics::ComputeStressTensor(nDim, tau_lam, gradvel, mu_lam);
  CNumerics::AddQCR(nDim, gradvel, tau_total, mu_turb / (mu_lam + mu_turb));
  CNumerics::AddQCR(nDim, gradvel, tau_turb, 1.0);
  for (size_t iDim = 0; iDim < nDim; iDim++)
    for (size_t jDim = 0; jDim < nDim; jDim++)
      REQUIRE(tau_total[iDim][jDim] == Approx(tau_lam[iDim][jDim] + tau_turb[iDim][jDim]).margin(1e-12));
}

namespace {

struct NEMOViscousFixture {
  static constexpr unsigned short nDim = 2;
  static constexpr unsigned short nSpecies = 5;
  static constexpr unsigned short nVar = nSpecies + nDim + 2;
  static constexpr unsigned short nPrimVar = nSpecies + nDim + 10;
  static constexpr unsigned short nPrimVarGrad = nSpecies + nDim + 8;
  static constexpr unsigned short T_INDEX = nSpecies;
  static constexpr unsigned short TVE_INDEX = nSpecies + 1;
  static constexpr unsigned short VEL_INDEX = nSpecies + 2;
  static constexpr unsigned short P_INDEX = nSpecies + nDim + 2;
  static constexpr unsigned short RHO_INDEX = nSpecies + nDim + 3;
  static constexpr unsigned short H_INDEX = nSpecies + nDim + 4;
  static constexpr unsigned short A_INDEX = nSpecies + nDim + 5;
  static constexpr unsigned short RHOCVTR_INDEX = nSpecies + nDim + 6;
  static constexpr unsigned short RHOCVVE_INDEX = nSpecies + nDim + 7;

  std::stringstream options;
  CConfig* config = nullptr;
  std::vector<su2double> primitive_i = std::vector<su2double>(nPrimVar, 0.0);
  std::vector<su2double> primitive_j = std::vector<su2double>(nPrimVar, 0.0);
  su2activematrix gradient_i = su2activematrix(nPrimVarGrad, nDim);
  su2activematrix gradient_j = su2activematrix(nPrimVarGrad, nDim);
  std::array<su2double, nSpecies> diffusion_i{};
  std::array<su2double, nSpecies> diffusion_j{};
  std::array<su2double, nSpecies> eve_i{};
  std::array<su2double, nSpecies> eve_j{};
  std::array<su2double, nSpecies> cvve_i{};
  std::array<su2double, nSpecies> cvve_j{};
  std::array<su2double, nVar> dT_i{};
  std::array<su2double, nVar> dT_j{};
  std::array<su2double, nVar> dTve_i{};
  std::array<su2double, nVar> dTve_j{};
  const std::array<su2double, nDim> coord_i{0.0, 0.0};
  const std::array<su2double, nDim> coord_j{3.0, 4.0};
  const std::array<su2double, nDim> normal{0.6, 0.8};

  NEMOViscousFixture() {
    options << "SOLVER= NEMO_NAVIER_STOKES\n"
            << "GAS_MODEL= AIR-5\n"
            << "GAS_COMPOSITION= (0.77, 0.23, 0.0, 0.0, 0.0)\n"
            << "FLUID_MODEL= SU2_NONEQ\n"
            << "FROZEN_MIXTURE= YES\n"
            << "TIME_DISCRE_FLOW= EULER_IMPLICIT\n"
            << "CONV_NUM_METHOD_FLOW= AUSM\n";
    config = new CConfig(options, SU2_COMPONENT::SU2_CFD, false);

    primitive_i[0] = primitive_j[0] = 0.77;
    primitive_i[1] = primitive_j[1] = 0.23;
    primitive_i[T_INDEX] = primitive_j[T_INDEX] = 300.0;
    primitive_i[TVE_INDEX] = primitive_j[TVE_INDEX] = 300.0;
    primitive_i[P_INDEX] = primitive_j[P_INDEX] = 101325.0;
    primitive_i[RHO_INDEX] = primitive_j[RHO_INDEX] = 1.0;
    primitive_i[H_INDEX] = primitive_j[H_INDEX] = 3.0e5;
    primitive_i[A_INDEX] = primitive_j[A_INDEX] = 340.0;
    primitive_i[RHOCVTR_INDEX] = primitive_j[RHOCVTR_INDEX] = 700.0;
    primitive_i[RHOCVVE_INDEX] = primitive_j[RHOCVVE_INDEX] = 1.0;
    gradient_i = su2double(0.0);
    gradient_j = su2double(0.0);
  }

  ~NEMOViscousFixture() { delete config; }

  void set_common(CNumerics& numerics) {
    numerics.SetCoord(coord_i.data(), coord_j.data());
    numerics.SetNormal(normal.data());
    numerics.SetPrimitive(primitive_i.data(), primitive_j.data());
    numerics.SetPrimVarGradient(CMatrixView<const su2double>(gradient_i), CMatrixView<const su2double>(gradient_j));
    numerics.SetDiffusionCoeff(diffusion_i.data(), diffusion_j.data());
    numerics.SetLaminarViscosity(2.0, 2.0);
    numerics.SetEddyViscosity(0.0, 0.0);
    numerics.SetThermalConductivity(0.0, 0.0);
    numerics.SetThermalConductivity_ve(0.0, 0.0);
    numerics.SetEve(eve_i.data(), eve_j.data());
    numerics.SetCvve(cvve_i.data(), cvve_j.data());
    numerics.SetdTdU(dT_i.data(), dT_j.data());
    numerics.SetdTvedU(dTve_i.data(), dTve_j.data());
  }

  void set_edge_linear_velocity_gradient() {
    gradient_i = su2double(0.0);
    gradient_j = su2double(0.0);
    constexpr su2double distance_squared = 25.0;
    for (unsigned short component = 0; component < nDim; ++component) {
      const su2double jump = primitive_j[VEL_INDEX + component] - primitive_i[VEL_INDEX + component];
      for (unsigned short dimension = 0; dimension < nDim; ++dimension) {
        const su2double edge_component = coord_j[dimension] - coord_i[dimension];
        gradient_i(VEL_INDEX + component, dimension) = jump * edge_component / distance_squared;
        gradient_j(VEL_INDEX + component, dimension) = jump * edge_component / distance_squared;
      }
    }
  }

  template <class Numerics>
  su2double directional_flux(Numerics& numerics) {
    set_common(numerics);
    const auto residual = numerics.ComputeResidual(config);
    su2double projected = 0.0;
    for (unsigned short component = 0; component < nDim; ++component)
      projected += normal[component] * residual.residual[nSpecies + component];
    return projected;
  }

  static su2double directional_jacobian(const su2double* const* matrix, const std::array<su2double, nDim>& direction) {
    su2double value = 0.0;
    for (unsigned short row = 0; row < nDim; ++row)
      for (unsigned short column = 0; column < nDim; ++column)
        value += direction[row] * matrix[nSpecies + row][nSpecies + column] * direction[column];
    return value;
  }
};

}  // namespace

TEST_CASE("NEMO corrected viscous residual returns distinct i and j Jacobians", "[NEMO][viscous][Jacobian]") {
  NEMOViscousFixture fixture;
  CAvgGradCorrected_NEMO numerics(fixture.nDim, fixture.nVar, fixture.nPrimVar, fixture.nPrimVarGrad, fixture.config);
  fixture.set_common(numerics);
  const auto base = numerics.ComputeResidual(fixture.config);
  const su2double analytic_i = fixture.directional_jacobian(base.jacobian_i, fixture.normal);
  const su2double analytic_j = fixture.directional_jacobian(base.jacobian_j, fixture.normal);

  REQUIRE(analytic_i == Approx(-8.0 / 15.0).epsilon(1.0e-12));
  REQUIRE(analytic_j == Approx(8.0 / 15.0).epsilon(1.0e-12));
  REQUIRE(analytic_i != Approx(analytic_j));
}
