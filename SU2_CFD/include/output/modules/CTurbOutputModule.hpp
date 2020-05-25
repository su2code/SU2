#pragma once

#include "COutputModule.hpp"

class CTurbOutputModule final : public CSolverOutputModule {

  int turb_model = NONE;
  int hybrid_RANSLES = NONE;
  int transModel = NONE;

public:
  explicit CTurbOutputModule(CConfig *config, int nDim) : CSolverOutputModule(nDim, config->GetKind_Turb_Model()),
    turb_model(config->GetKind_Turb_Model()),
    hybrid_RANSLES(config->GetKind_HybridRANSLES()),
    transModel(config->GetKind_Trans_Model()){}

  void LoadHistoryData(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                       const IterationInfo& iterationInfo) override;

  void DefineHistoryFields(CHistoryOutFieldManager& historyFields) override;

  void LoadVolumeData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                      const IterationInfo& iterationInfo, const PointInfo& pointInfo) override;

  void DefineVolumeFields(CVolumeOutFieldManager& volumeFields) override;
};

