#ifndef CFLOWCOEFFICIENTMODULE_HPP
#define CFLOWCOEFFICIENTMODULE_HPP

#include "COutputModule.hpp"

class CFlowCoefficientModule : public CSolverOutputModule {

public:
  explicit CFlowCoefficientModule(CConfig *config, int nDim);

  void DefineVolumeFields(CVolumeOutFieldManager& volumeFields) override;

  void LoadSurfaceData(CVolumeOutFieldManager& volumeFields) override;

};

#endif // CFLOWCOEFFICIENTMODULE_HPP
