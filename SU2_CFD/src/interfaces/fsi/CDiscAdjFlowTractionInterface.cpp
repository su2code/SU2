/*!
 * \file CDiscAdjFlowTractionInterface.cpp
 * \brief Declaration and inlines of the class to transfer flow tractions
 *        from a fluid zone into a structural zone in a discrete adjoint simulation.
 * \author Ruben Sanchez
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/interfaces/fsi/CDiscAdjFlowTractionInterface.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"
#include "../../../include/solvers/CSolver.hpp"

CDiscAdjFlowTractionInterface::CDiscAdjFlowTractionInterface(unsigned short val_nVar, unsigned short val_nConst,
                                                             const CConfig *config, bool conservative_) :
  CFlowTractionInterface(val_nVar, val_nConst, config, conservative_) {
}

void CDiscAdjFlowTractionInterface::GetPhysical_Constants(CSolver *flow_solution, CSolver *struct_solution,
                                                          CGeometry *flow_geometry, CGeometry *struct_geometry,
                                                          const CConfig *flow_config, const CConfig *struct_config){

  Preprocess(flow_config, struct_config, struct_geometry, struct_solution);

  if (!conservative) ComputeVertexAreas(struct_config, struct_geometry, struct_solution);

  /*--- No ramp applied ---*/
  Physical_Constants[1] = 1.0;
}
