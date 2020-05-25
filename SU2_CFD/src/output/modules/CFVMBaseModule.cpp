#include "../../../include/output/modules/CFVMBaseModule.hpp"

#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../include/solvers/CSolver.hpp"

void CFVMBaseModule::DefineVolumeFields(CVolumeOutFieldManager& volumeFields) {
  volumeFields.AddField("COORD_X", "x", "COORDINATES", "x-component of the coordinate vector", FieldType::DEFAULT);
  volumeFields.AddField("COORD_Y", "y", "COORDINATES", "y-component of the coordinate vector", FieldType::DEFAULT);
  if (nDim == 3)
    volumeFields.AddField("COORD_Z", "z", "COORDINATES", "z-component of the coordinate vector", FieldType::DEFAULT);

  volumeFields.AddField("NORMAL_X", "n_x", "NORMALS", "x-component of the normal vector", FieldType::DEFAULT);
  volumeFields.AddField("NORMAL_Y", "n_y", "NORMALS", "y-component of the normal vector", FieldType::DEFAULT);
  if (nDim == 3)
    volumeFields.AddField("NORMAL_Z", "n_z", "NORMALS", "z-component of the normal vector", FieldType::DEFAULT);

  volumeFields.AddField("AREA", "Area", "AREA", "Area", FieldType::SURFACE_INTEGRATE);

  if (gridMovement) {
    volumeFields.AddField("GRID_VELOCITY_X", "Grid_Velocity_x", "GRID_VELOCITY",
                          "x-component of the grid velocity vector", FieldType::DEFAULT);
    volumeFields.AddField("GRID_VELOCITY_Y", "Grid_Velocity_y", "GRID_VELOCITY",
                          "y-component of the grid velocity vector", FieldType::DEFAULT);
    if (nDim == 3)
      volumeFields.AddField("GRID_VELOCITY_Z", "Grid_Velocity_z", "GRID_VELOCITY",
                            "z-component of the grid velocity vector", FieldType::DEFAULT);
  }
}

void CFVMBaseModule::LoadSurfaceData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                                     const IterationInfo& iterationInfo, const PointInfo& pointInfo) {
  const auto* geometry = get<1>(solverData);
  const auto iVertex   = get<1>(pointInfo);
  const auto iMarker   = get<2>(pointInfo);
  const auto* normal = geometry->vertex[iMarker][iVertex]->GetNormal();

  su2double area = 0.0;
  for (int iDim = 0; iDim < nDim; iDim++) {
    area += normal[iDim] * normal[iDim];
  }
  area = sqrt(area);

  volumeFields.SetFieldValue("NORMAL_X", normal[0] / area);
  volumeFields.SetFieldValue("NORMAL_Y", normal[1] / area);
  if (nDim == 3) {
    volumeFields.SetFieldValue("NORMAL_Z", normal[2] / area);
  }
  volumeFields.SetFieldValue("AREA", area);
}

void CFVMBaseModule::LoadVolumeData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                                    const IterationInfo& iterationInfo, const PointInfo& pointInfo) {
  const auto* geometry  = get<1>(solverData);
  const auto  iPoint    = get<0>(pointInfo);
  const auto* coord   = geometry->nodes->GetCoord(iPoint);
  const auto* nodeGeo = geometry->nodes;

  volumeFields.SetFieldValue("COORD_X", coord[0]);
  volumeFields.SetFieldValue("COORD_Y", coord[1]);
  if (nDim == 3) volumeFields.SetFieldValue("COORD_Z", coord[2]);

  if (gridMovement) {
    volumeFields.SetFieldValue("GRID_VELOCITY_X", nodeGeo->GetGridVel()[iPoint][0]);
    volumeFields.SetFieldValue("GRID_VELOCITY_Y", nodeGeo->GetGridVel()[iPoint][1]);
    if (nDim == 3) volumeFields.SetFieldValue("GRID_VELOCITY_Z", nodeGeo->GetGridVel()[iPoint][2]);
  }
}