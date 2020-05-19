#pragma once

#include "COutputModule.hpp"

class CTurbOutputModule final : public CSolverOutputModule {

  int turb_model;

public:
  explicit CTurbOutputModule(CConfig *config, int nDim) : CSolverOutputModule(nDim, config->GetKind_Turb_Model()),
    turb_model(config->GetKind_Turb_Model()) {}

  void LoadHistoryData(CHistoryOutFieldManager& historyFields) override;

  void DefineHistoryFields(CHistoryOutFieldManager& historyFields) override;

  void LoadVolumeData(CVolumeOutFieldManager& volumeFields) override;

  void DefineVolumeFields(CVolumeOutFieldManager& volumeFields) override;
};

