/*!
 * \file ausm_tke_tests.cpp
 * \brief SLAU, SLAU2 and AUSM with the turbulent kinetic energy of SST in the total energy.
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
#include <array>
#include <cmath>
#include <memory>
#include <sstream>
#include "../../../SU2_CFD/include/numerics/flow/convection/ausm_slau.hpp"
#include "../../../SU2_CFD/include/numerics/flow/convection/fvs.hpp"

namespace {

constexpr su2double gamma = 1.4;

std::unique_ptr<CConfig> MakeConfig(bool implicit = false) {
  std::stringstream options;
  options << "SOLVER= EULER\nCONV_NUM_METHOD_FLOW= SLAU\nMACH_NUMBER= 0.5\n";
  options << "TIME_DISCRE_FLOW= "
          << (implicit ? "EULER_IMPLICIT\nUSE_ACCURATE_FLUX_JACOBIANS= YES\n" : "EULER_EXPLICIT\n");
  return std::make_unique<CConfig>(options, SU2_COMPONENT::SU2_CFD, false);
}

/*--- Primitive variables (T, u, v, [w], p, rho, h, c) of an ideal gas whose total energy contains k. ---*/
void Primitives(unsigned short nDim, su2double rho, const su2double* vel, su2double p, su2double k, su2double* V) {
  su2double q2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    V[iDim + 1] = vel[iDim];
    q2 += vel[iDim] * vel[iDim];
  }
  V[0] = p / (rho * 287.058);
  V[nDim + 1] = p;
  V[nDim + 2] = rho;
  V[nDim + 3] = gamma / (gamma - 1) * p / rho + 0.5 * q2 + k;
  V[nDim + 4] = sqrt(gamma * p / rho);
}

/*--- k changes neither the pressure nor the speed of sound: with the same density, velocity and pressure, the mass
 and momentum fluxes must not depend on k, and the energy flux changes only by the advected k. ---*/
void CheckTkeIndependence(CNumerics& numerics, const CConfig* config, unsigned short nDim, bool uniformTke = false) {
  const su2double vel_i[3] = {0.8, 0.1, -0.05}, vel_j[3] = {0.6, -0.1, 0.05}, n[3] = {2.0 / 7, 3.0 / 7, 6.0 / 7};
  const su2double n2[2] = {0.6, 0.8};
  su2double Vi[10] = {0.0}, Vj[10] = {0.0};
  numerics.SetNormal(nDim == 2 ? n2 : n);

  su2double flux[2][5] = {{0.0}};
  /*--- MSW mixes the two states, so k must be uniform for the fluxes to be independent of it. ---*/
  const su2double k_i[2] = {0.0, 0.3}, k_j[2] = {0.0, uniformTke ? 0.3 : 0.2};
  for (int c = 0; c < 2; c++) {
    Primitives(nDim, 1.2, vel_i, 1.1, k_i[c], Vi);
    Primitives(nDim, 0.9, vel_j, 0.8, k_j[c], Vj);
    numerics.SetPrimitive(Vi, Vj);
    numerics.SetTurbKineticEnergy(k_i[c], k_j[c]);
    const auto res = numerics.ComputeResidual(config);
    for (unsigned short iVar = 0; iVar < nDim + 2; iVar++) flux[c][iVar] = res.residual[iVar];
  }
  CAPTURE(nDim);
  for (unsigned short iVar = 0; iVar < nDim + 1; iVar++) CHECK(flux[1][iVar] == Approx(flux[0][iVar]).epsilon(1e-12));
  /*--- Energy: the extra flux is the mass flux times the upwind k (the flow is from i to j here). ---*/
  CHECK(flux[1][nDim + 1] - flux[0][nDim + 1] == Approx(flux[0][0] * k_i[1]).epsilon(1e-10));
}

/*--- Largest difference between the (accurate) Jacobians and central differences of the flux with respect to the
 conservative variables, k held fixed, relative to the largest entry. ---*/
su2double JacobianError(CNumerics& numerics, const CConfig* config, su2double k) {
  constexpr unsigned short nDim = 2, nVar = 4;
  using State = std::array<su2double, nVar>;
  const su2double normal[2] = {0.6, 0.8};
  numerics.SetNormal(normal);
  auto evaluate = [&](const State& Ui, const State& Uj, su2double* flux, su2double(*jac_i)[nVar],
                      su2double(*jac_j)[nVar]) {
    su2double Vi[10] = {0.0}, Vj[10] = {0.0};
    for (int side = 0; side < 2; side++) {
      const State& U = side ? Uj : Ui;
      const su2double vel[2] = {U[1] / U[0], U[2] / U[0]};
      const su2double p = (gamma - 1) * (U[3] - 0.5 * (U[1] * U[1] + U[2] * U[2]) / U[0] - U[0] * k);
      Primitives(nDim, U[0], vel, p, k, side ? Vj : Vi);
    }
    numerics.SetPrimitive(Vi, Vj);
    numerics.SetTurbKineticEnergy(k, k);
    const auto res = numerics.ComputeResidual(config);
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      flux[iVar] = res.residual[iVar];
      if (jac_i)
        for (unsigned short jVar = 0; jVar < nVar; jVar++) {
          jac_i[iVar][jVar] = res.jacobian_i[iVar][jVar];
          jac_j[iVar][jVar] = res.jacobian_j[iVar][jVar];
        }
    }
  };
  const State Ui = {1.2, 0.72, 0.24, 1.1 / (gamma - 1) + 0.5 * (0.72 * 0.72 + 0.24 * 0.24) / 1.2 + 1.2 * k};
  const State Uj = {0.9, 0.36, -0.18, 0.8 / (gamma - 1) + 0.5 * (0.36 * 0.36 + 0.18 * 0.18) / 0.9 + 0.9 * k};
  su2double flux[nVar], jac_i[nVar][nVar], jac_j[nVar][nVar];
  evaluate(Ui, Uj, flux, jac_i, jac_j);
  su2double maxErr = 0.0, maxJac = 0.0;
  for (int side = 0; side < 2; side++) {
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      State Up = side ? Uj : Ui, Um = Up;
      const su2double h = 1e-6;
      Up[jVar] += h;
      Um[jVar] -= h;
      su2double fp[nVar], fm[nVar];
      if (side) {
        evaluate(Ui, Up, fp, nullptr, nullptr);
        evaluate(Ui, Um, fm, nullptr, nullptr);
      } else {
        evaluate(Up, Uj, fp, nullptr, nullptr);
        evaluate(Um, Uj, fm, nullptr, nullptr);
      }
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        const su2double jac = side ? jac_j[iVar][jVar] : jac_i[iVar][jVar];
        maxErr = std::max(maxErr, std::abs(jac - (fp[iVar] - fm[iVar]) / (2 * h)));
        maxJac = std::max(maxJac, std::abs(jac));
      }
    }
  }
  return maxErr / maxJac;
}

}  // namespace

TEST_CASE("Accurate SLAU and AUSM+up Jacobians hold k fixed", "[AUSM][SST]") {
  auto config = MakeConfig(true);
  for (const su2double k : {0.0, 0.2}) {
    CAPTURE(k);
    CUpwSLAU_Flow slau(2, 4, config.get(), false);
    CHECK(JacobianError(slau, config.get(), k) < 1e-4);
    CUpwAUSMPLUSUP_Flow ausmup(2, 4, config.get());
    CHECK(JacobianError(ausmup, config.get(), k) < 1e-4);
  }
}

TEST_CASE("SLAU, SLAU2, AUSM, AUSM+up, AUSM+up2 and MSW do not use k in the speed of sound", "[AUSM][SST]") {
  auto config = MakeConfig();
  for (const unsigned short nDim : {2, 3}) {
    CUpwSLAU_Flow slau(nDim, nDim + 2, config.get(), false);
    CheckTkeIndependence(slau, config.get(), nDim);
    CUpwSLAU2_Flow slau2(nDim, nDim + 2, config.get(), false);
    CheckTkeIndependence(slau2, config.get(), nDim);
    CUpwAUSM_Flow ausm(nDim, nDim + 2, config.get());
    CheckTkeIndependence(ausm, config.get(), nDim);
    CUpwAUSMPLUSUP_Flow ausmup(nDim, nDim + 2, config.get());
    CheckTkeIndependence(ausmup, config.get(), nDim);
    CUpwAUSMPLUSUP2_Flow ausmup2(nDim, nDim + 2, config.get());
    CheckTkeIndependence(ausmup2, config.get(), nDim);
    CUpwMSW_Flow msw(nDim, nDim + 2, config.get());
    CheckTkeIndependence(msw, config.get(), nDim, true);
  }
}
