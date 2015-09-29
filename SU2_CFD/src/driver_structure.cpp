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
inline void CDriver::Run(CIteration **iteration_container,
                         COutput *output,
                         CIntegration ***integration_container,
                         CGeometry ***geometry_container,
                         CSolver ****solver_container,
                         CNumerics *****numerics_container,
                         CConfig **config_container,
                         CSurfaceMovement **surface_movement,
                         CVolumetricMovement **grid_movement,
                         CFreeFormDefBox*** FFDBox) { }


void CDriver::Preprocess(CIteration **iteration_container,
    COutput *output,
    CIntegration ***integration_container,
    CGeometry ***geometry_container,
    CSolver ****solver_container,
    CNumerics *****numerics_container,
    CConfig **config_container,
    CSurfaceMovement **surface_movement,
    CVolumetricMovement **grid_movement,
    CFreeFormDefBox*** FFDBox){

  unsigned short iZone;

  for (iZone = 0; iZone < nZone; iZone++) {
      Solver_Preprocessing(solver_container[iZone], geometry_container[iZone],
                               config_container[iZone], iZone);
      Integration_Preprocessing(integration_container[iZone], geometry_container[iZone],
                                    config_container[iZone], iZone);
      Numerics_Preprocessing(numerics_container[iZone], solver_container[iZone],
                                 geometry_container[iZone], config_container[iZone], iZone);
    }

}

CSingleZoneDriver::CSingleZoneDriver(CConfig **config, unsigned short val_nZone) : CDriver(config, val_nZone) { }

CSingleZoneDriver::~CSingleZoneDriver(void) { }

void CSingleZoneDriver::Run(CIteration **iteration_container,
                            COutput *output,
                            CIntegration ***integration_container,
                            CGeometry ***geometry_container,
                            CSolver ****solver_container,
                            CNumerics *****numerics_container,
                            CConfig **config_container,
                            CSurfaceMovement **surface_movement,
                            CVolumetricMovement **grid_movement,
                            CFreeFormDefBox*** FFDBox) {
  
  /*--- Run an iteration of the physics within this single zone.
   We assume that the zone of interest is in the ZONE_0 container position. ---*/
  
  iteration_container[ZONE_0]->Preprocess(); /*--- Does nothing for now. ---*/
  
  iteration_container[ZONE_0]->Iterate(output, integration_container, geometry_container,
                                       solver_container, numerics_container, config_container,
                                       surface_movement, grid_movement, FFDBox);
  
  iteration_container[ZONE_0]->Update(); /*--- Does nothing for now. ---*/
  
  iteration_container[ZONE_0]->Monitor(); /*--- Does nothing for now. ---*/
  
  iteration_container[ZONE_0]->Output(); /*--- Does nothing for now. ---*/
  
  iteration_container[ZONE_0]->Postprocess(); /*--- Does nothing for now. ---*/
  
}




CMultiZoneDriver::CMultiZoneDriver(CConfig **config, unsigned short val_nZone) : CDriver(config, val_nZone) {}

CMultiZoneDriver::~CMultiZoneDriver(void) { }

void CMultiZoneDriver::Run(CIteration **iteration_container,
                           COutput *output,
                           CIntegration ***integration_container,
                           CGeometry ***geometry_container,
                           CSolver ****solver_container,
                           CNumerics *****numerics_container,
                           CConfig **config_container,
                           CSurfaceMovement **surface_movement,
                           CVolumetricMovement **grid_movement,
                           CFreeFormDefBox*** FFDBox) {
  
  unsigned short iZone;
  
  /*--- Run a single iteration of a multi-zone problem by looping over all
   zones and executing the iterations. Note that data transers between zones
   and other intermediate procedures may be required. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    
    iteration_container[iZone]->Preprocess();  /*--- Does nothing for now. ---*/
    
    iteration_container[iZone]->Iterate(output, integration_container, geometry_container,
                                         solver_container, numerics_container, config_container,
                                         surface_movement, grid_movement, FFDBox);
    
    iteration_container[iZone]->Update();      /*--- Does nothing for now. ---*/
    
    iteration_container[iZone]->Monitor();     /*--- Does nothing for now. ---*/
    
    iteration_container[iZone]->Output();      /*--- Does nothing for now. ---*/
    
    iteration_container[iZone]->Postprocess(); /*--- Does nothing for now. ---*/
    
    //! TDE: here is where a call to CTransfer could be placed to move data
    //! between zones. Or, something similar for mixing planes, sliding meshes,
    //! time spectral, etc.
    
  }
  
}


CFSIDriver::CFSIDriver(CConfig **config, unsigned short val_nZone) : CDriver(config, val_nZone) {}

CFSIDriver::~CFSIDriver(void) { }

void CFSIDriver::Run(CIteration **iteration_container,
                     COutput *output,
                     CIntegration ***integration_container,
                     CGeometry ***geometry_container,
                     CSolver ****solver_container,
                     CNumerics *****numerics_container,
                     CConfig **config_container,
                     CSurfaceMovement **surface_movement,
                     CVolumetricMovement **grid_movement,
                     CFreeFormDefBox*** FFDBox) {
  
  //! TDE: For example, FSI BGS implementation could go here here.
  
}


