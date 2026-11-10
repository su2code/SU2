/*!
 * \file thermal_group_flux_tests.cpp
 * \brief Tests for thermal residuals at grouped boundary vertices.
 * \version 8.6.0 "Harrier"
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
 * Distributed under the GNU Lesser General Public License, version 2.1 or later.
 */
#include "catch.hpp"
#include "../../../SU2_CFD/include/solvers/thermal_group_flux.hpp"

TEST_CASE("Grouped wall thermal residual derivatives", "[BoundaryFaceGroups]") {
  const double area[] = {0.4, 0.7, 0.2}, distance[] = {0.1, 0.3, 0.5}, temperature[] = {1.1, 0.9, 0.8};
  auto group = [&](size_t i, double& a, double& d, double& t) { a = area[i]; d = distance[i]; t = temperature[i]; };
  auto flux = [&](double wall, double& derivative) {
    double value = 0.0;
    REQUIRE(ThermalGroupFlux::Isothermal(true, 3, 0.6, wall, group, value, derivative));
    return value;
  };
  const double wall = 1.2, step = 1e-6;
  double derivative = 0.0, scratch = 0.0;
  const double value = flux(wall, derivative);
  const double isoPlus = flux(wall + step, scratch);
  const double isoMinus = flux(wall - step, scratch);
  CHECK(derivative == Approx((isoPlus - isoMinus) / (2 * step)).epsilon(1e-8));
  double unused = 42.0;
  CHECK_FALSE(ThermalGroupFlux::Isothermal(false, 3, 0.6, wall, group, unused, scratch));
  CHECK(unused == 42.0);
  auto invalid = [&](size_t i, double& a, double& d, double& t) { group(i, a, d, t); if (i == 1) d = 0.0; };
  CHECK_FALSE(ThermalGroupFlux::Isothermal(true, 3, 0.6, wall, invalid, unused, scratch));
  auto transfer = [&](double t, double& jacobian) {
    double q = 0.0;
    ThermalGroupFlux::HeatTransfer(2.5, 1.7, t, 0.3, q, jacobian);
    return -q;
  };
  double transferJacobian = 0.0;
  const double transferResidual = transfer(wall, transferJacobian);
  const double transferPlus = transfer(wall + step, scratch);
  const double transferMinus = transfer(wall - step, scratch);
  CHECK(transferJacobian == Approx((transferPlus - transferMinus) / (2 * step)).epsilon(1e-8));
  const double rho = 1.3, momentum = 0.4, energy = 3.0, inverseCv = 0.4;
  auto conservativeResidual = [&](double density, double mom, double totalEnergy) {
    const double temperature = inverseCv * (totalEnergy / density - 0.5 * mom * mom / (density * density));
    double coefficient = 0.0;
    return transfer(temperature, coefficient);
  };
  const double temperatureWall = inverseCv * (energy / rho - 0.5 * momentum * momentum / (rho * rho));
  const double dTdrho = (-temperatureWall + inverseCv * 0.5 * momentum * momentum / (rho * rho)) / rho;
  const double densityFD = (conservativeResidual(rho + step, momentum, energy) -
                            conservativeResidual(rho - step, momentum, energy)) / (2 * step);
  const double energyFD = (conservativeResidual(rho, momentum, energy + step) -
                           conservativeResidual(rho, momentum, energy - step)) / (2 * step);
  CHECK(transferJacobian * dTdrho == Approx(densityFD).epsilon(1e-8));
  CHECK(transferJacobian * inverseCv / rho == Approx(energyFD).epsilon(1e-8));
  CHECK(derivative + transferJacobian == Approx((isoPlus + transferPlus - isoMinus - transferMinus) /
                                                 (2 * step)).epsilon(1e-8));
  CHECK(value + transferResidual == Approx(0.6 * ((wall - 1.1) * 4.0 + (wall - 0.9) * (0.7 / 0.3) +
                                                   (wall - 0.8) * 0.4) - 2.5 * (1.7 - wall) * 0.3));
}
