#pragma once
#include "COutputModule.hpp"

class CFVMBaseModule : public CSolverOutputModule {

  bool gridMovement = false;

public:
  explicit CFVMBaseModule(CConfig *config, int nDim) : CSolverOutputModule(nDim),
  gridMovement(config->GetGrid_Movement()){};

  void DefineVolumeFields(CVolumeOutFieldManager& volumeFields) override;

  void LoadSurfaceData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                       const IterationInfo& iterationInfo, const PointInfo& pointInfo) override;

  void LoadVolumeData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                      const IterationInfo& iterationInfo, const PointInfo& pointInfo) override;

};
