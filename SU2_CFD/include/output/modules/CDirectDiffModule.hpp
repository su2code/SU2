#pragma once

#include "COutputModule.hpp"

class CDirectDiffModule : public CSolverOutputModule {


  COutFieldCollection::InsertionVector modifiedFields, baseFields;

public:
  explicit CDirectDiffModule(CConfig *config);

  void DefineHistoryFieldModifier(COutFieldCollection& fieldCollection) override;

  void LoadHistoryDataModifier(COutFieldCollection& fieldCollection) override;

};
