/*!
 * \file roe_contact_tests.cpp
 * \brief Roe-type schemes must resolve a contact discontinuity and a shear layer exactly (upwind flux), in 2D
 *        and 3D, on faces not aligned with the axes.
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

std::unique_ptr<CConfig> MakeConfig(bool implicit) {
  std::stringstream options;
  options << "SOLVER= EULER\nCONV_NUM_METHOD_FLOW= ROE\nENTROPY_FIX_COEFF= 0.0\nROE_KAPPA= 0.5\n";
  options << "TIME_DISCRE_FLOW= " << (implicit ? "EULER_IMPLICIT" : "EULER_EXPLICIT") << "\n";
  return std::make_unique<CConfig>(options, SU2_COMPONENT::SU2_CFD, false);
}

/*--- Primitive (T, u, v, [w], p, rho, h, c) and secondary (dp/drho_e, dp/de_rho) variables of an ideal gas. ---*/
void State(unsigned short nDim, su2double rho, const su2double* vel, su2double p, su2double* V, su2double* S) {
  su2double q2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    V[iDim + 1] = vel[iDim];
    q2 += vel[iDim] * vel[iDim];
  }
  V[0] = p / (rho * 287.058);
  V[nDim + 1] = p;
  V[nDim + 2] = rho;
  V[nDim + 3] = gamma / (gamma - 1) * p / rho + 0.5 * q2;
  V[nDim + 4] = sqrt(gamma * p / rho);
  S[0] = p / rho;
  S[1] = (gamma - 1) * rho;
}

template <class Scheme>
std::unique_ptr<CNumerics> MakeScheme(unsigned short nDim, const CConfig* config) {
  return std::make_unique<Scheme>(nDim, nDim + 2, config);
}
template <>
std::unique_ptr<CNumerics> MakeScheme<CUpwRoe_Flow>(unsigned short nDim, const CConfig* config) {
  return std::make_unique<CUpwRoe_Flow>(nDim, nDim + 2, config, false);
}

/*--- Two states with the same pressure and normal velocity, different density and (with shear) tangential
 velocity: a contact discontinuity, whose exact (upwind) flux is the flux of the upstream state. ---*/
template <class Scheme>
void CheckContact(const char* name, unsigned short nDim, const su2double* vel, bool implicit, bool shear) {
  auto config = MakeConfig(implicit);
  auto numerics = MakeScheme<Scheme>(nDim, config.get());

  const su2double p = 1.0, rho_i = 1.4, rho_j = 0.7;
  const su2double n2[2] = {0.6, 0.8}, n3[3] = {2.0 / 7, 3.0 / 7, 6.0 / 7};
  const su2double* n = (nDim == 2) ? n2 : n3;
  const su2double t2[2] = {0.16, -0.12}, t3[3] = {0.15, -0.1, 0.0};  // tangential velocity jumps
  const su2double* t = (nDim == 2) ? t2 : t3;
  su2double vel_j[3] = {0.0};
  for (unsigned short iDim = 0; iDim < nDim; iDim++) vel_j[iDim] = vel[iDim] + (shear ? t[iDim] : 0.0);

  su2double Vi[10] = {0.0}, Vj[10] = {0.0}, Si[2], Sj[2];
  State(nDim, rho_i, vel, p, Vi, Si);
  State(nDim, rho_j, vel_j, p, Vj, Sj);
  numerics->SetPrimitive(Vi, Vj);
  numerics->SetSecondary(Si, Sj);
  numerics->SetNormal(n);
  const auto res = numerics->ComputeResidual(config.get());

  su2double un = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) un += vel[iDim] * n[iDim];
  const su2double rho_up = (un >= 0.0) ? rho_i : rho_j;
  const su2double* V_up = (un >= 0.0) ? Vi : Vj;

  CAPTURE(name, nDim, un, implicit, shear);
  CHECK(res.residual[0] == Approx(rho_up * un).margin(1e-12));
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    CHECK(res.residual[iDim + 1] == Approx(rho_up * V_up[iDim + 1] * un + p * n[iDim]).margin(1e-12));
  CHECK(res.residual[nDim + 1] == Approx(rho_up * V_up[nDim + 3] * un).margin(1e-12));
}

template <class Scheme>
void CheckContacts(const char* name, bool implicit = true, bool shear = true) {
  const su2double zero[3] = {0.0, 0.0, 0.0}, moving[3] = {0.3, -0.1, 0.2}, reverse[3] = {-0.3, 0.1, -0.2};
  for (const unsigned short nDim : {2, 3}) {
    for (const auto* vel : {zero, moving, reverse}) {
      CheckContact<Scheme>(name, nDim, vel, implicit, false);
      if (shear) CheckContact<Scheme>(name, nDim, vel, implicit, true);
    }
  }
}

}  // namespace

TEST_CASE("Roe schemes resolve a contact and a shear layer exactly", "[Roe]") {
  CheckContacts<CUpwRoe_Flow>("Roe");
  CheckContacts<CUpwGeneralRoe_Flow>("general Roe", true);
  CheckContacts<CUpwGeneralRoe_Flow>("general Roe", false);
  CheckContacts<CUpwLMRoe_Flow>("LMRoe");
  /*--- L2Roe scales the tangential velocity jump at low Mach number (not upwind by design). ---*/
  CheckContacts<CUpwL2Roe_Flow>("L2Roe", true, false);
}
