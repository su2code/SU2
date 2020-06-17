#pragma once
#include "COutputModule.hpp"

class CResidualModule : public CModifierModule {

  std::map<std::string, su2double> initialResiduals;
  const std::map<std::string, std::string> AverageGroupName = {{"BGS_RES", "bgs"},{"RMS_RES","rms"},{"MAX_RES", "max"}};
  bool TimeDomain;
public:
  explicit CResidualModule(CConfig *config, int nDim);

  void DefineHistoryFieldModifier(CHistoryOutFieldManager& historyFields) override;

  void LoadHistoryDataModifier(CHistoryOutFieldManager& historyFields,
                               const IterationInfo& iterationInfo) override;

};
