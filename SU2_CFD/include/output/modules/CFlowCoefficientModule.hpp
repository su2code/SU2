#ifndef CFLOWCOEFFICIENTMODULE_HPP
#define CFLOWCOEFFICIENTMODULE_HPP

#include "COutputModule.hpp"

class CFlowCoefficientModule : public CSolverOutputModule {

public:
  explicit CFlowCoefficientModule(CConfig *config);

  void DefineVolumeFields(COutFieldCollection& fieldCollection) override;

  void LoadSurfaceData(COutFieldCollection& fieldCollection) override;

};

#endif // CFLOWCOEFFICIENTMODULE_HPP
