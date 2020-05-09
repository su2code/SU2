#pragma once
#include "COutputModule.hpp"

class CResidualModule : public CSolverOutputModule {

  std::map<std::string, su2double> initialResiduals;
  const std::map<std::string, std::string> AverageGroupName = {{"BGS_RES", "bgs"},{"RMS_RES","rms"},{"MAX_RES", "max"}};

public:
  explicit CResidualModule(CConfig *config);

  void DefineHistoryFieldModifier(COutFieldCollection& fieldCollection) override;

  void LoadHistoryDataModifier(COutFieldCollection& fieldCollection) override;

};
