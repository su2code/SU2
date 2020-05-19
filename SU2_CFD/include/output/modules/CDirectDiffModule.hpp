#pragma once

#include "COutputModule.hpp"

class CDirectDiffModule : public CSolverOutputModule {


  COutFieldManager::FieldRefVector modifiedFields, baseFields;

public:
  explicit CDirectDiffModule(CConfig *config, int nDim);

  void DefineHistoryFieldModifier(CHistoryOutFieldManager& historyFields) override;

  void LoadHistoryDataModifier(CHistoryOutFieldManager& historyFields) override;

};
