/*!
 * \file output_baseline.cpp
 * \brief Main subroutines for flow discrete adjoint output
 * \author R. Sanchez
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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


#include "../../include/output/CBaselineOutput.hpp"

#include "../../../Common/include/geometry_structure.hpp"
#include "../../include/solver_structure.hpp"

CBaselineOutput::CBaselineOutput(CConfig *config, unsigned short nDim, CSolver* solver) : COutput(config, nDim, false) {
  
  unsigned short iField = 0;

  /*--- Set the requested volume fields to all fields in the solver ---*/
  
  requestedVolumeFields = solver->GetSolutionFields();
  
  /*--- remove the point ID from the fields --- */
  
  requestedVolumeFields.erase(requestedVolumeFields.begin());
  
  nRequestedVolumeFields = requestedVolumeFields.size();
    
  /*--- Remove first and last character of the strings (the quotation marks) ---*/
  
  for (iField = 0; iField < nRequestedVolumeFields; iField++){
    requestedVolumeFields[iField] = requestedVolumeFields[iField].substr(1, requestedVolumeFields[iField].size() - 2);
  }
     
  /*--- Set the volume filename --- */
  
  volumeFilename = "baseline";
  
  /*--- Set the surface filename ---*/
  
  surfaceFilename = "surface_baseline";
  
  /*--- Conv field --- */
  
  convField = "NONE";
  
}

CBaselineOutput::~CBaselineOutput(void) {


}

void CBaselineOutput::SetVolumeOutputFields(CConfig *config){

  unsigned short iField = 0;
    
  /*--- The first three fields should be the coordinates, if not, something is wrong ---*/
  
  if (requestedVolumeFields[0] != "x" || requestedVolumeFields[1] != "y"){
    SU2_MPI::Error("No coordinates found in the restart file!!", CURRENT_FUNCTION);
  }
  if (nDim == 3){
    if (requestedVolumeFields[2] != "z"){
      SU2_MPI::Error("No coordinates found in the restart file!!", CURRENT_FUNCTION);
    }
  }
  
  // Grid coordinates
  AddVolumeOutput(requestedVolumeFields[0], requestedVolumeFields[0], "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput(requestedVolumeFields[1], requestedVolumeFields[1], "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput(requestedVolumeFields[2], requestedVolumeFields[2], "COORDINATES", "z-component of the coordinate vector");

  // Add all the remaining fields
  
  for (iField = nDim; iField < nRequestedVolumeFields; iField++){
    AddVolumeOutput(requestedVolumeFields[iField], requestedVolumeFields[iField], "SOLUTION","");
  }
  
}

void CBaselineOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){
  
  unsigned short iField = 0;
  
  if (nRequestedVolumeFields != solver[0]->GetnVar()){
    SU2_MPI::Error("Number of requested fields and number of variables do not match.", CURRENT_FUNCTION);
  }
  
  /*--- Take the solver at index 0 --- */
  
  CVariable* Node_Sol  = solver[0]->node[iPoint];
   
  for (iField = 0; iField < nRequestedVolumeFields; iField++){
    SetVolumeOutputValue(requestedVolumeFields[iField], iPoint, Node_Sol->GetSolution(iField));
  }
  
}
