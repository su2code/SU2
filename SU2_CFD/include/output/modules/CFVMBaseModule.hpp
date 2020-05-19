#pragma once
#include "COutputModule.hpp"

class CFVMBaseModule : public CSolverOutputModule {


public:
  explicit CFVMBaseModule(CConfig *config, int nDim) : CSolverOutputModule(nDim){};

  void DefineVolumeFields(CVolumeOutFieldManager& volumeFields) override;

  void LoadSurfaceData(CVolumeOutFieldManager& volumeFields) override;

  void LoadVolumeData(CVolumeOutFieldManager& volumeFields) override;

};
