/*!
 * \file driver_structure.cpp
 * \brief The main subroutines for driving single or multi-zone problems.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 4.0.1 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../include/driver_structure.hpp"

CDriver::CDriver(CConfig **config, unsigned short val_nZone) {
  nZone = val_nZone;
}

CDriver::~CDriver(void) { }

//! TDE: Inline file for this perhaps? It's purely virtual.
inline void CDriver::Run(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                         CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                         CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) { }

CSingleZoneDriver::CSingleZoneDriver(CConfig **config, unsigned short val_nZone) : CDriver(config, val_nZone) { }

CSingleZoneDriver::~CSingleZoneDriver(void) { }

void CSingleZoneDriver::Run(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                            CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                            CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
  
  /*--- To be implemented after iteration overhaul. ---*/
  
}

CMultiZoneDriver::CMultiZoneDriver(CConfig **config, unsigned short val_nZone) : CDriver(config, val_nZone) {}

CMultiZoneDriver::~CMultiZoneDriver(void) { }

void CMultiZoneDriver::Run(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                           CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                           CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
  
  /*--- To be implemented after iteration overhaul. ---*/
  
}

CFSIDriver::CFSIDriver(CConfig **config, unsigned short val_nZone) : CDriver(config, val_nZone) {}

CFSIDriver::~CFSIDriver(void) { }

void CFSIDriver::Run(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                     CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
  
  /*--- For example, FSI BGS implementation here. ---*/
  
}

