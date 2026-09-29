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
#include <cmath>
#include <memory>
#include <sstream>
#include "../../../SU2_CFD/include/numerics/flow/convection/ausm_slau.hpp"

namespace {

constexpr su2double gamma = 1.4;

std::unique_ptr<CConfig> MakeConfig() {
  std::stringstream options;
  options << "SOLVER= EULER\nCONV_NUM_METHOD_FLOW= SLAU\nTIME_DISCRE_FLOW= EULER_EXPLICIT\n";
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
void CheckTkeIndependence(CNumerics& numerics, const CConfig* config, unsigned short nDim) {
  const su2double vel_i[3] = {0.8, 0.1, -0.05}, vel_j[3] = {0.6, -0.1, 0.05}, n[3] = {2.0 / 7, 3.0 / 7, 6.0 / 7};
  const su2double n2[2] = {0.6, 0.8};
  su2double Vi[10] = {0.0}, Vj[10] = {0.0};
  numerics.SetNormal(nDim == 2 ? n2 : n);

  su2double flux[2][5] = {{0.0}};
  const su2double k_i[2] = {0.0, 0.3}, k_j[2] = {0.0, 0.2};
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

}  // namespace

TEST_CASE("SLAU, SLAU2 and AUSM do not use k in the speed of sound", "[AUSM][SST]") {
  auto config = MakeConfig();
  for (const unsigned short nDim : {2, 3}) {
    CUpwSLAU_Flow slau(nDim, nDim + 2, config.get(), false);
    CheckTkeIndependence(slau, config.get(), nDim);
    CUpwSLAU2_Flow slau2(nDim, nDim + 2, config.get(), false);
    CheckTkeIndependence(slau2, config.get(), nDim);
    CUpwAUSM_Flow ausm(nDim, nDim + 2, config.get());
    CheckTkeIndependence(ausm, config.get(), nDim);
  }
}
