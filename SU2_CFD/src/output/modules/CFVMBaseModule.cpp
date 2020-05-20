#include "../../../include/output/modules/CFVMBaseModule.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../include/solvers/CSolver.hpp"

void CFVMBaseModule::DefineVolumeFields(CVolumeOutFieldManager& volumeFields){

  volumeFields.AddField("COORD_X", "x", "COORDINATES", "x-component of the coordinate vector", FieldType::DEFAULT);
  volumeFields.AddField("COORD_Y", "y", "COORDINATES", "y-component of the coordinate vector", FieldType::DEFAULT);
  if (nDim == 3)
    volumeFields.AddField("COORD_Z", "z", "COORDINATES", "z-component of the coordinate vector", FieldType::DEFAULT);

  volumeFields.AddField("NORMAL_X", "n_x", "NORMALS", "x-component of the normal vector", FieldType::DEFAULT);
  volumeFields.AddField("NORMAL_Y", "n_y", "NORMALS", "y-component of the normal vector", FieldType::DEFAULT);
  if (nDim == 3)
    volumeFields.AddField("NORMAL_Z", "n_z", "NORMALS", "z-component of the normal vector", FieldType::DEFAULT);

}

void CFVMBaseModule::LoadSurfaceData(CVolumeOutFieldManager &volumeFields){

  const auto* normal = solverData.vertex->GetNormal();
  su2double area = 0.0;
  for (int iDim = 0 ; iDim < nDim; iDim++){
    area += normal[iDim]*normal[iDim];
  }
  area = sqrt(area);

  volumeFields.SetFieldValue("NORMAL_X", normal[0]/area);
  volumeFields.SetFieldValue("NORMAL_Y", normal[1]/area);
  if (nDim == 3){
    volumeFields.SetFieldValue("NORMAL_Z", normal[2]/area);
  }
}

void CFVMBaseModule::LoadVolumeData(CVolumeOutFieldManager &volumeFields){

  const auto* coord = solverData.geometry->nodes->GetCoord(solverData.iPoint);

  volumeFields.SetFieldValue("COORD_X", coord[0]);
  volumeFields.SetFieldValue("COORD_Y", coord[1]);
  if (nDim == 3)
    volumeFields.SetFieldValue("COORD_Z", coord[2]);

}