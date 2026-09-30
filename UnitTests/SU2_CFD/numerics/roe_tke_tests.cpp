/*!
 * \file roe_tke_tests.cpp
 * \brief Roe eigenvector matrices with the turbulent kinetic energy of SST in the total energy.
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
#include "../../../SU2_CFD/include/numerics/flow/convection/roe.hpp"

namespace {

constexpr su2double gamma = 1.4;

std::unique_ptr<CConfig> MakeConfig(bool implicit = true) {
  std::stringstream options;
  options << "SOLVER= EULER\nCONV_NUM_METHOD_FLOW= ROE\nENTROPY_FIX_COEFF= 0.0\nROE_KAPPA= 0.5\n";
  options << "TIME_DISCRE_FLOW= " << (implicit ? "EULER_IMPLICIT" : "EULER_EXPLICIT") << "\n";
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

/*--- Secondary variables (dp/drho_e, dp/de_rho) of the ideal gas, for the general-gas schemes. ---*/
void Secondary(su2double rho, su2double p, su2double* S) {
  S[0] = p / rho;
  S[1] = (gamma - 1) * rho;
}

/*--- Common constructor for the schemes of the family. ---*/
template <class Scheme>
std::unique_ptr<CNumerics> MakeScheme(unsigned short nDim, unsigned short nVar, const CConfig* config) {
  return std::make_unique<Scheme>(nDim, nVar, config);
}
template <>
std::unique_ptr<CNumerics> MakeScheme<CUpwRoe_Flow>(unsigned short nDim, unsigned short nVar, const CConfig* config) {
  return std::make_unique<CUpwRoe_Flow>(nDim, nVar, config, false);
}

/*--- Roe flux between two states with the same pressure and normal velocity and different density, k and (with
 shear) tangential velocity. With the entropy fix off this is a contact discontinuity, whose exact (upwind) flux is
 the flux of the upstream state. ---*/
template <class Scheme>
void CheckContactWithTkeJump(const char* name, unsigned short nDim, const su2double* vel, bool implicit, bool shear) {
  auto config = MakeConfig(implicit);
  const unsigned short nVar = nDim + 2;
  auto numerics = MakeScheme<Scheme>(nDim, nVar, config.get());

  const su2double p = 1.0, rho_i = 1.4, rho_j = 0.7, k_i = 0.1, k_j = 0.3;
  const su2double n2[2] = {0.6, 0.8}, n3[3] = {2.0 / 7, 3.0 / 7, 6.0 / 7};
  const su2double* n = (nDim == 2) ? n2 : n3;
  const su2double t2[2] = {0.16, -0.12}, t3[3] = {0.15, -0.1, 0.0};  // tangential velocity jumps
  const su2double* t = (nDim == 2) ? t2 : t3;
  su2double vel_j[3] = {0.0};
  for (unsigned short iDim = 0; iDim < nDim; iDim++) vel_j[iDim] = vel[iDim] + (shear ? t[iDim] : 0.0);
  su2double Vi[10] = {0.0}, Vj[10] = {0.0}, Si[2], Sj[2];
  Primitives(nDim, rho_i, vel, p, k_i, Vi);
  Primitives(nDim, rho_j, vel_j, p, k_j, Vj);
  Secondary(rho_i, p, Si);
  Secondary(rho_j, p, Sj);
  numerics->SetPrimitive(Vi, Vj);
  numerics->SetSecondary(Si, Sj);
  numerics->SetNormal(n);
  numerics->SetTurbKineticEnergy(k_i, k_j);
  const auto res = numerics->ComputeResidual(config.get());

  su2double un = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) un += vel[iDim] * n[iDim];
  const su2double rho_up = (un >= 0.0) ? rho_i : rho_j;
  const su2double* V_up = (un >= 0.0) ? Vi : Vj;

  CAPTURE(name, nDim, un, implicit, shear);
  CHECK(res.residual[0] == Approx(rho_up * un).margin(1e-12));
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    CHECK(res.residual[iDim + 1] == Approx(rho_up * V_up[iDim + 1] * un + p * n[iDim]).margin(1e-12));
  CHECK(res.residual[nVar - 1] == Approx(rho_up * V_up[nDim + 3] * un).margin(1e-12));
}

template <class Scheme>
void CheckContactsWithTkeJump(const char* name, bool implicit = true, bool shear = true) {
  const su2double zero[3] = {0.0, 0.0, 0.0}, moving[3] = {0.3, -0.1, 0.2}, reverse[3] = {-0.3, 0.1, -0.2};
  for (const unsigned short nDim : {2, 3}) {
    for (const auto* vel : {zero, moving, reverse}) {
      CheckContactWithTkeJump<Scheme>(name, nDim, vel, implicit, false);
      if (shear) CheckContactWithTkeJump<Scheme>(name, nDim, vel, implicit, true);
    }
  }
}

}  // namespace

TEST_CASE("Roe schemes resolve a contact with a jump of k exactly", "[Roe][SST]") {
  CheckContactsWithTkeJump<CUpwRoe_Flow>("Roe");
  CheckContactsWithTkeJump<CUpwGeneralRoe_Flow>("general Roe", true);
  CheckContactsWithTkeJump<CUpwGeneralRoe_Flow>("general Roe", false);
  CheckContactsWithTkeJump<CUpwLMRoe_Flow>("LMRoe");
  /*--- L2Roe scales the tangential velocity jump at low Mach number (not upwind by design). ---*/
  CheckContactsWithTkeJump<CUpwL2Roe_Flow>("L2Roe", true, false);
}

TEST_CASE("Roe eigenvector matrices are inverse with k in the total energy", "[Roe][SST]") {
  auto config = MakeConfig();
  for (const unsigned short nDim : {2, 3}) {
    const unsigned short nVar = nDim + 2;
    CUpwRoe_Flow numerics(nDim, nVar, config.get(), false);
    const su2double velocity[3] = {0.4, -0.3, 0.2}, normal[3] = {2.0 / 7, 3.0 / 7, 6.0 / 7};
    const su2double n2[3] = {0.6, 0.8, 0.0};
    su2double P[5][5], Pinv[5][5];
    numerics.GetPMatrix(1.1, velocity, 1.2, nDim == 3 ? normal : n2, P, 0.3);
    numerics.GetPMatrix_inv(1.1, velocity, 1.2, nDim == 3 ? normal : n2, Pinv, 0.3);
    for (unsigned short i = 0; i < nVar; i++) {
      for (unsigned short j = 0; j < nVar; j++) {
        su2double prod = 0.0;
        for (unsigned short l = 0; l < nVar; l++) prod += P[i][l] * Pinv[l][j];
        CAPTURE(nDim, i, j);
        CHECK(prod == Approx(i == j ? 1.0 : 0.0).margin(1e-12));
      }
    }
  }
}

TEST_CASE("Roe does not dissipate a stationary contact with k in the total energy", "[Roe][SST]") {
  /*--- Zero velocity, same pressure and k on both sides, different density: a contact discontinuity at rest.
   The Roe flux must be the pressure only (no mass or energy flux) when the entropy fix is off. ---*/
  auto config = MakeConfig();
  constexpr unsigned short nDim = 2, nVar = nDim + 2;
  CUpwRoe_Flow numerics(nDim, nVar, config.get(), false);

  const su2double zero[2] = {0.0, 0.0}, p = 1.0, k = 0.2;
  su2double Vi[10] = {0.0}, Vj[10] = {0.0}, normal[2] = {0.6, 0.8};
  Primitives(nDim, 1.4, zero, p, k, Vi);
  Primitives(nDim, 0.7, zero, p, k, Vj);
  numerics.SetPrimitive(Vi, Vj);
  numerics.SetNormal(normal);
  numerics.SetTurbKineticEnergy(k, k);
  const auto res = numerics.ComputeResidual(config.get());

  CHECK(res.residual[0] == Approx(0.0).margin(1e-12));
  CHECK(res.residual[1] == Approx(p * normal[0]).epsilon(1e-12));
  CHECK(res.residual[2] == Approx(p * normal[1]).epsilon(1e-12));
  CHECK(res.residual[3] == Approx(0.0).margin(1e-12));
}
