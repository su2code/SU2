#pragma once

#include "COutputModule.hpp"

class CTurbOutputModule final : public CSolverOutputModule {

  int turb_model;

public:
  explicit CTurbOutputModule(CConfig *config) : CSolverOutputModule(config->GetKind_Turb_Model()),
    turb_model(config->GetKind_Turb_Model()) {}

  void LoadHistoryData(COutFieldCollection& fieldCollection) override;

  void DefineHistoryFields(COutFieldCollection& fieldCollection) override;

};

