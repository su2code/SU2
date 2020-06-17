#pragma once

#include "COutputModule.hpp"

class CDirectDiffModule : public CModifierModule {


  COutFieldManager::FieldRefVector modifiedFields, baseFields;

public:
  explicit CDirectDiffModule(CConfig *config, int nDim);

  void DefineHistoryFieldModifier(CHistoryOutFieldManager& historyField) override;

  void LoadHistoryDataModifier(CHistoryOutFieldManager& historyFields,
                               const IterationInfo& iterationInfo) override;

};
