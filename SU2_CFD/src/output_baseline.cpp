#include "../include/output/output_baseline.hpp"

CBaselineOutput::CBaselineOutput(CConfig *config, CGeometry *geometry, CSolver *solver, unsigned short val_iZone) : COutput(config) {

  nDim = geometry->GetnDim();
  
  unsigned short iField = 0;

  /*--- Set the requested volume fields to all fields in the solver ---*/
  
  RequestedVolumeFields = solver->GetSolutionFields();
  
  /*--- remove the point ID from the fields --- */
  
  RequestedVolumeFields.erase(RequestedVolumeFields.begin());
  
  nRequestedVolumeFields = RequestedVolumeFields.size();
    
  /*--- Remove first and last character of the strings (the quotation marks) ---*/
  
  for (iField = 0; iField < nRequestedVolumeFields; iField++){
    RequestedVolumeFields[iField] = RequestedVolumeFields[iField].substr(1, RequestedVolumeFields[iField].size() - 2);
  }
     
  /*--- Set the volume filename --- */
  
  VolumeFilename = "baseline";
  
  /*--- Set the surface filename ---*/
  
  SurfaceFilename = "surface_baseline";
  
}

CBaselineOutput::~CBaselineOutput(void) {


}

void CBaselineOutput::SetVolumeOutputFields(CConfig *config){

  unsigned short iField = 0;
    
  /*--- The first three fields should be the coordinates, if not, something is wrong ---*/
  
  if (RequestedVolumeFields[0] != "x" || RequestedVolumeFields[1] != "y"){
    SU2_MPI::Error("No coordinates found in the restart file!!", CURRENT_FUNCTION);
  }
  if (nDim == 3){
    if (RequestedVolumeFields[2] != "z"){
      SU2_MPI::Error("No coordinates found in the restart file!!", CURRENT_FUNCTION);
    }
  }
  
  // Grid coordinates
  
  AddVolumeOutput(RequestedVolumeFields[0], RequestedVolumeFields[0], "COORDINATES");
  AddVolumeOutput(RequestedVolumeFields[1], RequestedVolumeFields[1], "COORDINATES");
  if (nDim == 3)
    AddVolumeOutput(RequestedVolumeFields[2], RequestedVolumeFields[2], "COORDINATES");

  // Add all the remaining fields
  
  for (iField = nDim-1; iField < nRequestedVolumeFields; iField++){
    AddVolumeOutput(RequestedVolumeFields[iField], RequestedVolumeFields[iField], "SOLUTION");
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
    SetVolumeOutputValue(RequestedVolumeFields[iField], iPoint, Node_Sol->GetSolution(iField));
  }
  
}