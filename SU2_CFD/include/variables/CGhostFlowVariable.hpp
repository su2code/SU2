/*!
 * \file CGhostFlowVariable.hpp
 * \brief Flow variables of the ghost points of one scalar solver's boundary marker.
 * \author P. Gomes
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

#include "CFlowVariable.hpp"

/*!
 * \class CGhostFlowVariable
 * \brief Flow variables of the ghost points of one marker, indexed by vertex.
 * \note Sized to the largest marker and to the primitive layout of the flow solver, so one set
 *       of indices reads the ghost points and the interior ones. Only the primitives the flux
 *       kernels read are filled by the boundary's fill pass (density, velocity, laminar and
 *       eddy viscosity); the gradients, the limiters and the non-physical edge counter are not,
 *       and the first order boundary path does not read them.
 */
class CGhostFlowVariable final : public CFlowVariable {
 public:
  CGhostFlowVariable(unsigned long npoint, unsigned long ndim, unsigned long nvar, unsigned long nprimvar,
                     unsigned long nprimvargrad, const CConfig* config)
      : CFlowVariable(npoint, ndim, nvar, nprimvar, nprimvargrad, config) {}

  /*!
   * \brief Never read: a ghost is not part of the dual-time residual, only its primitives are.
   */
  inline su2double GetDensity_time_n(unsigned long) const override { return 0.0; }
  inline su2double GetDensity_time_n1(unsigned long) const override { return 0.0; }
};
