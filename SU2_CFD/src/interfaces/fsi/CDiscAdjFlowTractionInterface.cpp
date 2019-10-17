/*!
 * \file CDiscAdjFlowTractionInterface.cpp
 * \brief Declaration and inlines of the class to transfer flow tractions
 *        from a fluid zone into a structural zone in a discrete adjoint simulation.
 * \author Ruben Sanchez
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

CDiscAdjFlowTractionInterface::CDiscAdjFlowTractionInterface(void) : CFlowTractionInterface() {

}

CDiscAdjFlowTractionInterface::CDiscAdjFlowTractionInterface(unsigned short val_nVar, unsigned short val_nConst,
                                                             CConfig *config) :
  CFlowTractionInterface(val_nVar, val_nConst, config) {

}

CDiscAdjFlowTractionInterface::~CDiscAdjFlowTractionInterface(void) {

}

void CDiscAdjFlowTractionInterface::GetPhysical_Constants(CSolver *flow_solution, CSolver *struct_solution,
                                                          CGeometry *flow_geometry, CGeometry *struct_geometry,
                                                          CConfig *flow_config, CConfig *struct_config){

  /*--- We have to clear the traction before applying it, because we are "adding" to node and not "setting" ---*/

  struct_solution->GetNodes()->Clear_FlowTraction();

  Preprocess(flow_config);

  /*--- No ramp applied ---*/
  Physical_Constants[1] = 1.0;
}
