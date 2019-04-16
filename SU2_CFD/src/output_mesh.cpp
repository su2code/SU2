

#include "../include/output/output_mesh.hpp"

CMeshOutput::CMeshOutput(CConfig *config, CGeometry *geometry, unsigned short val_iZone) : COutput(config) {

  nDim = geometry->GetnDim();

  /*--- Set the default history fields if nothing is set in the config file ---*/
  
  RequestedVolumeFields.push_back("COORDINATES");
  nRequestedVolumeFields = RequestedVolumeFields.size();
    
  /*--- Set the volume filename --- */
  
  VolumeFilename = config->GetMesh_Out_FileName();
  
  /*--- Set the surface filename ---*/
  
  SurfaceFilename = "surface_mesh";
  
}

CMeshOutput::~CMeshOutput(void) {


}

void CMeshOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES");

  
}

void CMeshOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CPoint*    Node_Geo  = geometry->node[iPoint];
 
  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));  
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));
  
}
